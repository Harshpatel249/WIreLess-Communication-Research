import csv
import numpy as np
import pandas as pd
from matplotlib import pyplot as plt


#Reading the data from csv files
pairpwr = pd.read_csv('/content/drive/My Drive/WCProject/data/1.csv')    #car_id, RSU_id, Received_pwr
RSUload = pd.read_csv('/content/drive/My Drive/WCProject/data/2.csv')    #Worst case RSUload
pairts = pd.read_csv('/content/drive/My Drive/WCProject/data/3.csv')     #No of links per time-step

#Converting the data to numpy 2-d array for further manipulation
pairpwr = np.asarray(pairpwr)
RSUload = np.asarray(RSUload)
pairts = np.asarray(pairts)

max = 0 #Upper limit of records for each time-step
min = 0 #Lower limit of records for each time-step
itr = 0 #Keep track of iterations
totalhandovers = 0 #Total handovers for the simulation

effRSULoad = np.zeros_like([(0,0,0)],shape=(int(len(RSUload)/len(pairts)),1)) #RSUload/pairts gives number of RSUs
finalopall = [] #To store output of all timesteps

#For storing effective load distribution
RSU1load = np.zeros_like([(0)],shape=(200,))
RSU2load = np.zeros_like([(0)],shape=(200,))
RSU3load = np.zeros_like([(0)],shape=(200,))

#For stroing worst case distribution
RSU1loadw = np.zeros_like([(0)],shape=(200,))
RSU2loadw = np.zeros_like([(0)],shape=(200,))
RSU3loadw = np.zeros_like([(0)],shape=(200,))

for i in pairts:        #iterating through no of pairs per time-step
    #Getting records for single time-step
    num = int(i)
    min = max
    max = max + num
    arr = np.array(pairpwr[min:max,:])   #Getting all the records for current time-step
    arr = arr[arr[:,0].argsort()]        #Sorting by car id   
    
    #Processing the records for given timep-step
    if(itr>0):
       finalopall.append(finalop)
    #Variables for output of current time-step
    effRSULoad[0] = 0
    effRSULoad[1] = 0
    effRSULoad[2] = 0
    finalop = [] #List to store final carid,RSUid and received power for current time-step
    pcar_id = -1 #previous car_id (for checking if the vehicle can be connected to multiple RSUs)
    pRSU = -1 #previous RSU (for changing the load distribution)
    ppwr = 1 #previous received power

    handovers = 0 #No of handovers
    print("itr:")
    print(itr)
  
    #Processing
    for k in arr:
        car_id = int(k[0])
        car_RSU = int(k[1])
        car_pwr = float(k[2])
        if(car_id != pcar_id):
            finalop.append([car_id,car_RSU,car_pwr])
            effRSULoad[(car_RSU-1)] = (int(effRSULoad[(car_RSU-1)])+1)  #Keeping track of load on RSU
            pcar_id = car_id   #Book keeping
            pRSU = car_RSU     #Book keeping
            ppwr = car_pwr     #book keeping
            
        elif(car_id == pcar_id):    #if car is connected to multiple RSU
            if(ppwr<-90 or (car_pwr-ppwr)>10):              #if previous link is not feasible anymore or new link has power differnce of 10Db
                if(car_pwr>-90):
                    finalop.pop()
                    finalop.append([car_id,car_RSU,car_pwr])
                    effRSULoad[(car_RSU-1)] = (int(effRSULoad[(car_RSU-1)])+1)   #Keeping track of load on RSU
                    effRSULoad[(pRSU-1)] = (int(effRSULoad[(pRSU-1)])-1)         #Keeping track of load on RSU
            elif((car_pwr>-90 and (car_pwr-ppwr)>0) and (effRSULoad[pRSU-1]>effRSULoad[car_RSU-1] and effRSULoad[pRSU-1]>10)):      #Equal RSU load distribution
                    finalop.pop()
                    finalop.append([car_id,car_RSU,car_pwr])
                    effRSULoad[(car_RSU-1)] = (int(effRSULoad[(car_RSU-1)])+1)   #Keeping track of load on RSU
                    effRSULoad[(pRSU-1)] = (int(effRSULoad[(pRSU-1)])-1)         #Keeping track of load on RSU
            
    
    print("RSU Load: \n")
    print(effRSULoad)
    RSU1load[itr] = effRSULoad[0]
    RSU2load[itr] = effRSULoad[1]
    RSU3load[itr] = effRSULoad[2]
    
    
    print("Car_id, RSU_id, Pwr \n")
    for p in finalop:
        print(p) 
    
    
    #Calculating handovers
    if(itr>0):
        for m in finalop:
            for n in finalopall[itr-1]:
                if(m[0] == n[0]):
                    if(m[1] != n[1]):
                        handovers = handovers + 1
                elif(n[0]>m[0]):
                    break
                
        print("Handover for current time-step: ")
        print(handovers)
        totalhandovers = totalhandovers + handovers
    
 
    itr = itr + 1    #book keeping

print("Total handovers: ")
print(totalhandovers)

#Plotting
'''
#Effective load distribution
print(RSU1load)
print(RSU2load)
print(RSU3load)
'''
c = 0 #book-keeping for RSU1
v = 0 #book-keeping for RSU1
b = 0 #book-keeping for RSU1

for i in RSUload:
    if(i[0]==1):
        RSU1loadw[c] = i[1]
        c = c + 1
    elif(i[0]==2):
        RSU2loadw[v] = i[1]
        v = v + 1
    elif(i[0]==3):
        RSU3loadw[b] = i[1]
        b = b + 1
'''
#Worst case load
print(RSU1loadw)
print(RSU2loadw)
print(RSU3loadw)
'''
xaxis = range(0,200)

plt.plot(xaxis,RSU1loadw)
plt.plot(xaxis,RSU1load)
plt.legend(["Worst-Case load","Effective load distribution"])
plt.xlabel("Time-steps (seconds)")
plt.ylabel("No of vehicles connected(load)")
plt.title("Load distribution for RSU:1")

plt.plot(xaxis,RSU2loadw)
plt.plot(xaxis,RSU2load)
plt.legend(["Worst-Case load","Effective load distribution"])
plt.xlabel("Time-steps (seconds)")
plt.ylabel("No of vehicles connected(load)")
plt.title("Load distribution for RSU:2")

plt.plot(xaxis,RSU3loadw)
plt.plot(xaxis,RSU3load)
plt.legend(["Worst-Case load","Effective load distribution"])
plt.xlabel("Time-steps (seconds)")
plt.ylabel("No of vehicles connected(load)")
plt.title("Load distribution for RSU:3")