__author__ = "filip"
__version__ = "2018.11.30"
                        
import rhinoscriptsyntax as rs
import random




BQ = len(BUILDINGS_NAMES) # BUILDINGS Quantity - the amount of buildings
BQ2 = BQ**2
BC = BUILDINGS_STARTINGCOORDINATES

BUILDINGS_Coordinates_inTIME = [] # lists of coordinates of all buildings for all iterations
BUILDINGS_CORES_Coordinates_inTIME = []
BUILDINGS_CORES_Colors_inTIME = []
BUILDINGS_Occupancy_inTIME = [] # lists of occupancy of all the buildings for all iterations

UQ = len(USERS) # USER Quantity - the amount of users
UC=[] # users corinates
USER_X_coordinates = []
USERS_Colors_inTIME=[] # colors for all users in all iterations
USERS_Coordinates_inTIME = []

PODS_Coordinates_inTIME = []
PODS_Colors_inTIMES = []
PODS_Occupancy_inTIME = []

PATHS_inTIME =[]
avrageCOMMUTEdistance_inTIME =[]
avrageCOMMUTEtime_inTIME =[]

users_ACTIVITIES_inTIME = []
TEXT_coordinate_inTIME = []
TEXT_inTIME = [] # to delete!
PRINT_LOG_inTIME = []
PRINT_LOG = [str(),str(),str(),str(),str()]
TIME = []

iNR = 0 
Vzero = rs.CreateVector(0,0,0)
Pzero = rs.CreatePoint(0,0,0)
hbS = boxSIZE/2 # half of box size
pheryphery_border1 = squeezing_area * hbS
pheryphery_border2 = boxSIZE - pheryphery_border1
squeezing_multiplier = squeezingPOWER/(pheryphery_border1)**2
maxHALOsize = BUILDINGS_Halo
maxflow = 0


def IDcomponent(cooridante):
    id_component, distance = divmod(cooridante,module)
    id_component = int(id_component)
    id_component = min(id_component, partitions-1)
    id_component = max(id_component, 0)
    return id_component


def Squeezer(coordinate):
    if coordinate < pheryphery_border1:
        squeezerVECTORcoordinate = min((coordinate - pheryphery_border1)**2*squeezing_multiplier,BUILDINGS_Halo)
    elif coordinate > pheryphery_border2:
        squeezerVECTORcoordinate = -min((coordinate - pheryphery_border2)**2*squeezing_multiplier,BUILDINGS_Halo)
    else:
        squeezerVECTORcoordinate = 0
    return squeezerVECTORcoordinate


##################################################################################
#################################### SKY CITY zoneS ##############################
##################################################################################
                        
module = boxSIZE / partitions
N = [-1,0,1]
ZONES = dict()
                        
class Zone(object):
                        
    def __init__(self, x,y,z):
        self.id = 'C%.2d%.2d%.2d' % (x, y, z)
        self.coordinates = module*(x+0.5), module*(y+0.5), module*(z+0.5)
        self.neighbours = []
        self.buildings = []
        self.users = []
        self.pods = []
    def __repr__(self):
        return self.id
                        
                        
for x in range(partitions):
    for y in range(partitions):
        for z in range(partitions):
            zone = Zone(x,y,z)
            key = zone.id
            ZONES[key] = zone
            for xn in N:
                for yn in N:
                    for zn in N:
                        xf = x + xn
                        yf = y + yn
                        zf = z + zn
                        if xf >= 0 and yf >=0 and zf >=0 and xf < partitions and yf < partitions and zf < partitions:
                            ZONES[key].neighbours = ZONES[key].neighbours+['C%.2d%.2d%.2d' % (xf, yf, zf)]


##################################################################################
#################################### CELLS #######################################
##################################################################################
                        
print '>>> max size of a building: ', max(BUILDINGS_SIZES)
size_of_a_grid_for_buildings = 3*max(BUILDINGS_SIZES)
core_sizes = [1]
                        
                        
GC_Plus = [False, False, False]
GC_Size = [0,0,0]
GC_Rows = [[0],[0],[0]]
CELLS = [rs.CreateVector(0,0,0)]
gridcells = 1
gvc = dict()
a = 0
                        
while gridcells < size_of_a_grid_for_buildings: 
    b = (a + 1) % 3
    c = (b + 1) % 3
    if GC_Plus[a]:
        GC_Size[a] = -GC_Size[a]
        GC_Plus[a] = False
        GC_Rows[a] = [GC_Size[a]]+GC_Rows[a]
    else:
        GC_Size[a] = -GC_Size[a]+1
        GC_Plus[a] = True
        GC_Rows[a] = GC_Rows[a] + [GC_Size[a]]
                            
    gvc['c'+str(a)] = [GC_Size[a]]
    gvc['c'+str(b)] = GC_Rows[b]
    gvc['c'+str(c)] = GC_Rows[c]
                            
    for x in gvc['c0']:
        for y in gvc['c1']:
            for z in gvc['c2']:
                CELLS.append(rs.CreateVector(x,y,z))
                gridcells += 1
    core_sizes.append(gridcells)
    a = b
                        
for c in range(len(CELLS)):
    CELLS[c] = CELLS[c]*4
                        
def Core_Size(amount_of_pods):
    s = int(amount_of_pods/5)+1
    counter = 0
    core_size = core_sizes[counter]
    while core_sizes[counter] < s:
        counter += 1
        core_size = core_sizes[counter]
    return core_size
                        
##################################################################################
################################  !BUILDINGS!   ##################################
##################################################################################
                        
class Building(object):
                        
    def __init__(self, id, point, color, size):
        self.id = id
        self.coordinates = point
        self.zone = 'C%.2d%.2d%.2d' % (IDcomponent(point[0]), IDcomponent(point[1]), IDcomponent(point[2]))
        self.color = color
        self.size = size
        self.members = []
        self.users = []
        self.pods = []
        self.flows = dict()
        self.flowsum = 0
        self.flowvector = Vzero
        self.separationvector = Vzero
        self.occupancy = 0
        self.movementvector = Vzero
        
        cs = Core_Size(self.size)
        self.cellsofcore = CELLS[0:cs]
        self.cubes = CELLS[cs:cs+size*3]

    def __zone__(self):
        x = IDcomponent(self.coordinates[0])
        y = IDcomponent(self.coordinates[1])
        z = IDcomponent(self.coordinates[2])
        self.zone = 'C%.2d%.2d%.2d' % (x, y, z)

    def __addflow__(self,buildig):
        self.flowsum += 1
        if buildig not in self.flows:
            self.flows[buildig] = 1
        else:
            self.flows[buildig] += 1

    def __reserve_cell__(self):
        cube = self.cubes[random.randint(0,len(self.cubes)-1)]
        self.cubes.remove(cube)
        return cube

    def  __release_from__(self,pod):
        self.cubes.append(pod.cell)

    def __move_to__(self,b):
        rout = ROUTS[self.id+b]
        self.flowvector = rs.VectorAdd(rout.unit * self.flows[b], self.flowvector)

    def __separate_from__(self,B):
        self.separationvector = self.separationvector + ROUTS[B+self.id].unit * min(BUILDINGS_Halo/2,BUILDINGS_Halo**2/ROUTS[B+self.id].distance**2)
        if abs(self.separationvector[0]) > 10 or abs(self.separationvector[0]) > 10 or abs(self.separationvector[0]) > 10:
            print self.separationvector

    def __clean__(self):
        self.__zone__()
        self.flowvector = Vzero
        self.separationvector = Vzero
        self.movementvector = Vzero

    def __totalclean__(self):
        self.occupancy = 0
        self.flows = dict()
        self.flowsum = 0
        global maxflow
        maxflow = 0
                        
### FOR BUILDING CLASS
##  FOR ZONE def
#        if x in bounderyindecies or y in bounderyindecies or z in bounderyindecies:
#            self.peripheral = True
#        else:
#            self.peripheral = False
##  FOR __addflow__ def
#        global maxflow 
#        if self.flowsum > maxflow:
#            maxflow = self.flowsum
                        
BQ = len(BUILDINGS_NAMES)
BUILDINGS = dict()
BUILDINGS_List = []
for b in range(BQ):
    building = Building(BUILDINGS_NAMES[b],BUILDINGS_STARTINGCOORDINATES[b],BUILDINGS_COLORS[b],BUILDINGS_SIZES[b])
    BUILDINGS[building.id] = building
    BUILDINGS_List.append(building)
    for c in ZONES[building.zone].neighbours:
        ZONES[c].buildings.append(building.id)


##################################################################################
#################################### !PODS! ######################################
##################################################################################
                        
activities = ['HOUSE ACTIVITY','SLEEP','WORK','SOCIALIZING','SPORTS','SHOPPING','CIVIC','EDUCATION','RELIGION']
ACTIVITIES_COLORS = dict()
for a in range(9):
    ACTIVITIES_COLORS[activities[a]] = activities_colors[a]

PODS_dict = dict()
avaliblePODS = dict()
for a in activities:
    avaliblePODS[a] =[]

PQ = 0 # pods quanitity


class Pod(object):

    def __init__(self,PQ,user):
        self.id = 'P{0:05}'.format(PQ)
        self.activity = user.nextactivity
        self.color = ACTIVITIES_COLORS[self.activity]
        self.user = user.nr
        self.nextlocation = user.nextlocation
        self.departed = False
        self.todiscard = False
        x = random.randint(0,boxSIZE)
        y = random.randint(0,boxSIZE)
        height = boxSIZE*2
        self.coordinates = rs.CreatePoint(x,y,height)

    def __reserve__(self,cell):
        self.nextcell = cell

    def __update_plan__(self):
        self.departed = False
        self.todiscard = False
        self.location = self.nextlocation
        self.cell = self.nextcell
        return user.nr[1:0]

    def __discard_after_that_slot__(self):
        self.todiscard = True

    def __fly__(self,progress):
        self.destination = rs.PointAdd(BUILDINGS[self.nextlocation].coordinates,self.nextcell)
        vector = rs.VectorCreate(self.destination,self.coordinates)/leftiterations
        self.coordinates = rs.PointAdd(self.coordinates,vector)

    def __stay__(self):
        self.coordinates = rs.PointAdd(BUILDINGS[self.location].coordinates,self.cell)

    def __sink__(self,modifier):
        height = -(self.coordinates[2]-1.8)
        v = rs.CreateVector(0,0,height/boxSIZE*2*BUILDINGS_Speed/modifier)
        self.coordinates = rs.PointAdd(self.coordinates,v)

PODS_list = []

for user in USERS:
    user.spatium = 0.0
    user.iterations_of_commuting = 0
    user.activity = user.schedule_activities[starting_schedule_slot]
    user.nextactivity = user.activity
    user.location = user.schedule_locations[starting_schedule_slot]
    user.nextlocation = user.location

    pod = Pod(PQ,user)
    key = pod.id
    PODS_list.append(pod)
    PODS_dict[key] = pod
    PQ +=1


    building = BUILDINGS[user.nextlocation]
    cell = building.__reserve_cell__()
    pod.__reserve__(cell)
    pod.coordinates = rs.PointAdd(building.coordinates,pod.nextcell)
    user.coordinates = pod.coordinates
    user.__reserve_pod__(pod)
    user.pod = user.nextpod



#####################################################################################
#####################################  !ROUTS!   ####################################
#####################################################################################
                        
class Route(object):
                        
    def __init__(self, b1, b2):
        self.id = b1+b2
        self.counted = False
        self.vector = Vzero
        self.distance = 0.0
        self.unit = Vzero
                        
    def __update__(self, other,b1,b2):
                                
        self.counted = True
        other.counted = True
                                
        self.vector = rs.VectorCreate(BUILDINGS[b2].coordinates,BUILDINGS[b1].coordinates)
        other.vector = - self.vector
                                
        self.distance = rs.VectorLength(self.vector)
        other.distance = self.distance
                                
        self.unit = self.vector / self.distance
        other.unit = - self.unit
                        
    def __clean__(self):
        self.counted = False
                        
    def __repr__(self):
        return self.id


ROUTS = dict()
for b1 in BUILDINGS:
    for b2 in BUILDINGS:
        if b1 != b2:
            routename = b1+b2
            ROUTS[routename] = Route(b1,b2)
                        
print 'amount of routs:', len(ROUTS)
                        
                        
##################################################################################
'################################################################################'
################################# SIMMULATION ####################################
'################################################################################'
'################################################################################'
################################# !START! ########################################
'################################################################################'
################## SCHEDUAL predition & redtrodiction Slots ######################
'################################################################################'

iterationduration = int(900/iterations_per_SCHEDUAL_Slot)
schedule_lenght = len(USERS[0].schedule_locations)
amount_of_iterations = simmulation_lenght * iterations_per_SCHEDUAL_Slot
simulations_slots = range(starting_schedule_slot, starting_schedule_slot+simmulation_lenght)
tested_slots = range(starting_schedule_slot+retrodiction, starting_schedule_slot+simmulation_lenght+prediction)
prediction_and_retrodicion_slots = range(retrodiction,prediction+1)
amount_of_transitions = len(prediction_and_retrodicion_slots)-1
                        
for s in range(simmulation_lenght):
    simulations_slots[s] = simulations_slots[s] % schedule_lenght
for s in range(simmulation_lenght+prediction-retrodiction):
    if tested_slots > 0:
        tested_slots[s] = tested_slots[s] % schedule_lenght


for frame in range(simmulation_lenght):

    slot = simulations_slots[frame]
    nextslot = (slot + 1) % schedule_lenght
    slotrange = []
    for s in prediction_and_retrodicion_slots:
        slotrange.append(tested_slots[frame-retrodiction+s])

    for pod in PODS_list:
        if pod.todiscard:
            BUILDINGS[pod.nextlocation].__release_from__(pod)
            avaliblePODS[pod.activity].append(pod)
        pod.__update_plan__()

    Mutable_Users = []

    for user in USERS:
        user.changes_activity = False
        for t in range(amount_of_transitions):
            if user.schedule_locations[slotrange[t]] != user.schedule_locations[slotrange[t+1]]:
                BUILDINGS[user.schedule_locations[slotrange[t]]].__addflow__(user.schedule_locations[slotrange[t+1]])
                BUILDINGS[user.schedule_locations[slotrange[t+1]]].__addflow__(user.schedule_locations[slotrange[t]])
        user.location = user.schedule_locations[slot]
        user.activity = user.schedule_activities[slot]
        user.nextlocation = user.schedule_locations[nextslot]
        user.nextactivity = user.schedule_activities[nextslot]
        user.color = ACTIVITIES_COLORS[user.activity]
        user.pod = user.nextpod
        previous_pod = user.pod
        previous_pod.nextlocation = user.nextlocation

        if user.activity != user.nextactivity:
            Mutable_Users.append(user)
            user.changes_activity = True
            building = BUILDINGS[user.nextlocation]
            nextcell = building.__reserve_cell__()

            if user.location != user.nextlocation:
                user.relocating = True
                user.departed = False

                user.relocating = True
                user.departed = False
                user.podtravel = random.choice([True,False])
                if user.podtravel:
                    temporarycell = building.__reserve_cell__()
                    previous_pod.__reserve__(temporarycell)
                    previous_pod.__discard_after_that_slot__()
                    user.text = '{0} travels by pod to {1} object!'.format(user.first_name,user.nextactivity)
                else:
                    user.text = '{0} goes to {1} object!'.format(user.first_name,user.nextactivity)
            else:
                user.relocating = False
                previous_pod.__discard_after_that_slot__()
                user.text = '{0} asks for {1} pod!'.format(user.first_name,user.nextactivity)

            if avaliblePODS[user.nextactivity] == []:
                users_pod = Pod(PQ,user)
                key = users_pod.id
                PODS_dict[key] = users_pod
                PODS_list.append(users_pod)
                PQ += 1
            else:
                emptypodsTYPE = avaliblePODS[user.nextactivity]
                users_pod = random.choice(emptypodsTYPE)
                emptypodsTYPE.remove(users_pod)

            users_pod.__reserve__(nextcell)
            users_pod.nextlocation = user.nextlocation
            user.nextpod = users_pod
        else:
            user.relocating = False
            


    print 'amount of pods:', PQ

    print_log_counter = 0
    mutable_users_amount = len(Mutable_Users)
    text_frequency = mutable_users_amount / iterations_per_SCHEDUAL_Slot
    print_moments = []
    if mutable_users_amount > 0:
        for i in range(mutable_users_amount):
            print_moments.append(i*text_frequency)


#################################################################################################
########################################## ITERATIONS ###########################################
#################################################################################################


#    print("Simmulation slot: "+str(frame)+"    |    schedual slots number: "+str(slot))
    for iteration in range(iterations_per_SCHEDUAL_Slot):

        s_count = int(iNR*iterationduration)
        m_count, seconds = divmod(s_count, 60)
        hour_counter, minute_counter = divmod(m_count, 60)
        nychthemeron, hour_counter = divmod(hour_counter, 24)
        nychthemeron += 1

        minutes = m_count + starting_schedule_slot*15
        hours, minutes = divmod(minutes, 60)
        day, hours = divmod(hours, 24)
#        if nychthemeron == 1:
        time = '{0:01}:{1:02}:{2:02}'.format(hours, minutes, seconds)
#        else:
#            time = '{0} nychthemeron {1:01}:{2:02}:{3:02}'.format(nychthemeron, hours, minutes, seconds)

        print 'TIME:', time,"    |    Simmulation slot: "+str(frame)+"    |    schedual Slots number: "+str(slot)+"    |    iteration of a slot: "+str(iteration)+"    |    iteration: "+str(iNR)


        BUILDINGS_squeezingVECTORS = [] # separation vectors between all buildings 
        BUILDINGS_MOVEMNT=[]


        for b in BUILDINGS:
            ZONES[BUILDINGS[b].zone].buildings = []
        for b in BUILDINGS:
            for c in ZONES[BUILDINGS[b].zone].neighbours:
                ZONES[c].buildings.append(BUILDINGS[b].id)
        for r in ROUTS:
            ROUTS[r].__clean__()


######################## VECTORS and DISTANCES between all OBJECTS ########################
                        
                        
        for b1 in BUILDINGS:
            b = BUILDINGS[b1]
            for b2 in BUILDINGS[b1].flows:
                if ROUTS[b1+b2].counted == False:
                    ROUTS[b1+b2].__update__(ROUTS[b2+b1], b1, b2)
                b.__move_to__(b2)
            b.flowvector = b.flowvector / max(b.flowsum,1)
                                    
                        
            for b2 in ZONES[b.zone].buildings:
                if b2 != b1:
                    if ROUTS[b1+b2].counted == False:
                        ROUTS[b1+b2].__update__(ROUTS[b2+b1], b1, b2)
                    b.__separate_from__(b2)

                        
                        
        for b in BUILDINGS_NAMES:
            bc = BUILDINGS[b].coordinates
#            if BUILDINGS[BUILDINGS_NAMES[i]].peripheral == True:
            BUILDINGS_squeezingVECTORS.append(rs.CreateVector(Squeezer(bc[0]),Squeezer(bc[1]),Squeezer(bc[2]-5)))
#            else:
#                BUILDINGS_squeezingVECTORS.append(Vzero)
                        
                        
######################## MOVEMNET of OBJECTS according to DISTANCES ########################

        BUILDINGS_CORES_Coordinates_inTIME.append([])
        BUILDINGS_CORES_Colors_inTIME.append([])
        BUILDINGS_Coordinates_inTIME.append([])
        
        modifier = max(int(iterations_per_SCHEDUAL_Slot/10)-iteration,1)
        
        for i in range(BQ):
            b = BUILDINGS[BUILDINGS_NAMES[i]]
            BUILDINGS_MOVEMNT.append(Vzero)
            BUILDINGS_MOVEMNT[i]=rs.VectorAdd(BUILDINGS_MOVEMNT[i],BUILDINGS_squeezingVECTORS[i])
            BUILDINGS_MOVEMNT[i]=rs.VectorAdd(BUILDINGS_MOVEMNT[i],b.flowvector)
            BUILDINGS_MOVEMNT[i]=rs.VectorAdd(BUILDINGS_MOVEMNT[i],b.separationvector)
            b.movementvector = BUILDINGS_MOVEMNT[i]*BUILDINGS_Speed/modifier 
                
            b.coordinates = rs.VectorAdd(BC[i],b.movementvector)
            BC[i] = b.coordinates
            BUILDINGS_Coordinates_inTIME[iNR].append(b.coordinates)
            for cell in b.cellsofcore:
                BUILDINGS_CORES_Coordinates_inTIME[iNR].append(rs.VectorAdd(b.coordinates,cell))
                BUILDINGS_CORES_Colors_inTIME[iNR].append(b.color)


################################# MOVEMNET of USERS #################################



        leftiterations = float(iterations_per_SCHEDUAL_Slot-iteration)
        USERS_Coordinates_inTIME.append([])
        USERS_Colors_inTIME.append([])

        if iteration < 2.0/3.0*iterations_per_SCHEDUAL_Slot:
            for user in USERS:
                if user.activity != user.nextactivity:
                    newpod = user.nextpod
                    newpod.__fly__(leftiterations)
                pod_in_use = user.pod
                pod_in_use.__stay__()
                user.coordinates = pod_in_use.coordinates
                USERS_Coordinates_inTIME[iNR].append(user.coordinates)
                USERS_Colors_inTIME[iNR].append(user.color) #C

                
        else:
            migrationprogress = round(3.0*float(iteration)/float(iterations_per_SCHEDUAL_Slot)-2,2)
            print 'iNR:', iNR, 'slot simulation progress:', migrationprogress,  '%'
            for user in USERS:
                departurelocation = BUILDINGS[user.location]
                oldpod = user.pod
                newpod = user.nextpod
                if user.relocating:
                    if oldpod.departed == False:
                        vector = rs.VectorCreate(BUILDINGS[user.nextlocation].coordinates,user.coordinates)

                        if migrationprogress > 2.0/3.0*(1.0-rs.VectorLength(vector)/dis):
                            newpod.__fly__(leftiterations)
                            oldpod.departed = True
                            departurelocation.__release_from__(oldpod)
                            if user.podtravel == True:
                                oldpod.__fly__(leftiterations)
                                user.__follow__(oldpod,leftiterations)
                                print 'user', user.first_name, user.last_name, 'left his object'
                            else:
                                oldpod.__stay__()
                                avaliblePODS[oldpod.activity].append(oldpod) 
                                user.__go_to__(newpod.destination,leftiterations)
                        else:
                            newpod.__fly__(leftiterations)
                            oldpod.__stay__()
                            user.coordinates = oldpod.coordinates

                    else:
                        newpod.__fly__(leftiterations)
                        if user.podtravel:
                            oldpod.__fly__(leftiterations)
                            user.__follow__(oldpod,leftiterations)
                        else:
                            user.__go_to__(newpod.destination,leftiterations)
                else:
                    if user.activity != user.nextactivity:
                        newpod.__fly__(leftiterations)
                    oldpod.__stay__()
                    user.coordinates = oldpod.coordinates
                
                USERS_Coordinates_inTIME[iNR].append(user.coordinates)
                USERS_Colors_inTIME[iNR].append(user.color) #C
        
        for podtype in avaliblePODS:
            if avaliblePODS[podtype] != []:
                for emptypod in avaliblePODS[podtype]:
                    emptypod.__sink__(modifier)



        PODS_Coordinates_inTIME.append([])
        PODS_Colors_inTIMES.append([])
        for pod in PODS_list:
            PODS_Coordinates_inTIME[iNR].append(pod.coordinates)
            PODS_Colors_inTIMES[iNR].append(pod.color)


        commute_distance = 0.0
        total_iterations_of_commuting = 0
        PATHS_inTIME.append(([],[]))
        for user in USERS:
            commute_distance += user.spatium
            total_iterations_of_commuting += user.iterations_of_commuting
            if user.relocating:
                PATHS_inTIME[iNR][0].append(user.coordinates)
                PATHS_inTIME[iNR][1].append(BUILDINGS[user.nextlocation].coordinates)

        avarage_commute_time_in_seconds = int(total_iterations_of_commuting/UQ*iterationduration)
        commute_min, commute_sec = divmod(avarage_commute_time_in_seconds,60)
        commute_h, commute_min = divmod(commute_min,60)
        if commute_h > 0:
            avarage_commute_time ='{0} h {1} m {2} s'.format(commute_h, commute_min, commute_sec)
        elif commute_min > 0:
            avarage_commute_time ='{0} min {1} sec'.format(commute_min, commute_sec)
        else:
            avarage_commute_time ='{0} seconds'.format(commute_sec)

        avarage_commute_distance_in_meters = round(commute_distance/UQ,2)
        commute_km, commute_m = divmod(int(avarage_commute_distance_in_meters),1000)
        if commute_km > 0:
            avarage_commute_distance ='{0} km {1} m'.format(commute_km, commute_m)
        else:
            avarage_commute_distance ='{0} meters'.format(avarage_commute_distance_in_meters)
        avrageCOMMUTEdistance_inTIME.append(avarage_commute_distance)
        avrageCOMMUTEtime_inTIME.append(avarage_commute_time)
        TIME.append(time)

        activities_amounts = []
        for a in range(9):
            activity = activities[a]
            activities_amounts.append(0)
            for u in USERS:
                if u.activity == activity:
                    activities_amounts[a] +=1
        users_ACTIVITIES_inTIME.append(activities_amounts)
        PODS_Occupancy_inTIME.append(round(UQ/PQ,2))


        if mutable_users_amount > 0:
            while iteration >= print_moments[print_log_counter] and print_log_counter < mutable_users_amount:
                PRINT_LOG = PRINT_LOG[-5:]
                PRINT_LOG.append(Mutable_Users[print_log_counter].text)
                print_log_counter +=1
 
        TEXT_inTIME.append([])
        TEXT_coordinate_inTIME.append([])
        for user in USERS:
            if user.changes_activity == True: 
                TEXT_inTIME[iNR].append(user.text)
                TEXT_coordinate_inTIME[iNR].append(user.coordinates)

        iNR +=1
        for b in BUILDINGS:    
            BUILDINGS[b].__clean__()
    for b in BUILDINGS:    
        BUILDINGS[b].__totalclean__()

USERS_commute_data = []
for user in USERS:
    u_commute_sec = user.iterations_of_commuting*iterationduration
    u_commute_min, u_commute_sec = divmod(u_commute_sec,60)
    u_commute_h, u_commute_min = divmod(u_commute_min,60)
    USERS_commute_data.append([user.first_name,user.last_name,int(user.spatium),'{0:02}:{1:02}:{2:02}'.format(u_commute_h, u_commute_min, u_commute_sec)])