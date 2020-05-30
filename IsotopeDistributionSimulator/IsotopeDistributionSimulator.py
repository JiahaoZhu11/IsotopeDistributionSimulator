# Created By: Jiahao Zhu
# Created Date: 2020/5/30

import math
import matplotlib.pyplot as plt

print("============== Isotope Distribution Simulator ==============")

# function for showing the command list

def help():
    print("\n----------------------------------------\n")
    print('Input "exit" to quit')
    print('Input "edit" to edit isotope list')
    print('Input "help" to display the command list')

help()

# function for commands

def Operation(command, isotopes, editing = False):
    if command == "exit":
        exit(0)
    elif command == "edit":
        if editing:
            print("Editing!");
        else:
            editIso(isotopes)
    elif command == "help":
        help()
    else:
        return False
    return True
    

# function for adding isotopes

def editIso(isotopes):
    print("\n----------------------------------------\n")
    print('Please input isotopes in the format of "Name, Mass, Relative Intensity"')
    print('Input "file" to import data from file')
    print('Input "clean" to clean all isotopes')
    print('Input "finish" to terminate editing')
    while True:
        iso = input("Isotope " + str(len(isotopes) + 1) + ": ").replace(" ","")
        while Operation(iso, isotopes,True):
            print ("\n----------------------------------------\n")
            iso = input("Isotope " + str(len(isotopes) + 1) + ": ").replace(" ","")
        if iso == "finish":
            if isotopes:
                break
            print("Please add some isotopes before continue!")
        elif iso == "clean":
            isotopes = []
            continue
        elif iso == "file":
            path = input("File Path: ")
            try:
                with open(path) as lines:
                    for line in lines:
                        isotopes.append((line.split(",")[0].replace('"','').replace("'",""),
                                        float(line.split(",")[1]),float(line.split(",")[2])))
            except Exception:
                print("Invalid file!")
        else:
            try:
                isotopes.append((iso.split(",")[0].replace('"','').replace("'",""),
                                float(iso.split(",")[1]),float(iso.split(",")[2])))
            except Exception:
                print("Invalid input!")
    return

# initial a list of isotopes

isotopes = []
editIso(isotopes)

while True:

    # obtain compound name & make a dictionary of atoms included

    print("\n----------------------------------------\n")
    compName = input("Compound (e.g. C_1_H_1_Cl_3): ").replace(" ","")
    while Operation(compName, isotopes):
        print("\n----------------------------------------\n")
        compName = input("Compound (e.g. C_1_H_1_Cl_3): ").replace(" ","")
    compound = compName.split("_")

    compoundDict = {}
    atomName = ""

    for key in compound:
        if atomName == "":
            atomName = key
        else:
            try:
                compoundDict[atomName] = int(key)
            except Exception:
                print("Invalid compound name!")
                continue
            atomName = ""

    # function for calculating all the coefficient of the compound

    def compCal(isotopes, compName, compoundDict):
        print("Calculation Proceeding...")
        compProp = [compName.replace('_','')]
        atomProp = atomCal(isotopes,compoundDict)
        if not atomProp:
            return ()
        print("Calculating terms for the compound...")
        for atom in compoundDict:
            compCalHelper(atomProp, atom, compProp)
        return tuple(compProp)
    
    def compCalHelper(atomProp, atom, compProp):
        if len(compProp) == 1:
            return compProp.extend(list(atomProp[atom]))
        else:
            newTerms = {}
            for mass1 in compProp[1]:
                for mass2 in atomProp[atom][0]:
                    newTerms[mass1 + mass2] = compProp[1][mass1] * atomProp[atom][0][mass2]
            compProp[1] = newTerms
            compProp[2] = compProp[2] + atomProp[atom][1]
            compProp[3] = math.sqrt(math.pow(compProp[3], 2) + math.pow(atomProp[atom][2], 2))

    # function for calculating all the coefficient of each atom
    
    def atomCal(isotopes,compoundDict):
        atomProp = {}
        print("Calculating terms for each element...")
        for atom in compoundDict:
            print("Calculating for " + atom + "...")
            terms = {}
            if compoundDict[atom] > 0:
                selIso = []
                isoComb = []
                for iso in isotopes:
                    if iso[0][-len(atom):] == atom:
                        selIso.append(iso)
                        if not isoComb:
                            isoComb.append(1)
                        else:
                            isoComb.append(0)
                if not selIso:
                    print(atom + " is not included in the isotope list!")
                    return []
                indexDict = {tuple(isoComb): selIso[0][2]}
                isoComb[0] -= 1
                aveMass = selIso[0][1] * selIso[0][2]
                for i in range(1, len(isoComb)):
                    tempComb = list(isoComb)
                    tempComb[i] += 1
                    indexDict[tuple(tempComb)] = selIso[i][2]
                    aveMass += selIso[i][1] * selIso[i][2]
                SD = 0
                for iso in selIso:
                    SD += iso[2] * math.pow(iso[1] - aveMass, 2)
                aveMass = aveMass * compoundDict[atom]
                SD = math.sqrt(SD * compoundDict[atom])
                if compoundDict[atom] > 1:
                    indexDict = atomCalHelper(selIso, compoundDict[atom] - 1, indexDict)
                for index in indexDict:
                    mass = 0
                    for i in range(len(index)):
                        mass += index[i] * selIso[i][1]
                    terms[mass] = indexDict[index]
                atomProp[atom] = (terms, aveMass, SD)
        return atomProp

    def atomCalHelper(selIso, numOfAtom, indexDict):
        newIndexDict = {}
        if numOfAtom == 1:
            for isoComb in indexDict:
                for i in range (len(selIso)):
                    tempComb = list(isoComb)
                    tempComb[i] += 1
                    if tuple(tempComb) in newIndexDict:
                        newIndexDict[tuple(tempComb)] += indexDict[isoComb] * selIso[i][2]
                    else:
                        newIndexDict[tuple(tempComb)] = indexDict[isoComb] * selIso[i][2]
            return newIndexDict
        if numOfAtom > 1:
            for isoComb in indexDict:
                for i in range (len(selIso)):
                    tempComb = list(isoComb)
                    tempComb[i] += 1
                    if tuple(tempComb) in newIndexDict:
                        newIndexDict[tuple(tempComb)] += indexDict[isoComb] * selIso[i][2]
                    else:
                        newIndexDict[tuple(tempComb)] = indexDict[isoComb] * selIso[i][2]
            return atomCalHelper(selIso, numOfAtom - 1, newIndexDict) 

    # obtain data for the compound

    compProp = compCal(isotopes, compName, compoundDict)
    if not compProp:
        continue
    print("Collecting data...")
    terms = compProp[1]
    aveMass = compProp[2]
    SD = compProp[3]
    massRange = 10 * math.sqrt(1 + math.pow(SD,2))

    massList = []
    intensityList = []
    maxIntensity =  0
    totalIntensity = 0

    for mass in terms:
        massList.append(mass)
        intensityList.append(terms[mass])
        if maxIntensity < terms[mass]:
            maxIntensity = terms[mass]
        totalIntensity += terms[mass]

    if abs(totalIntensity - 1) > 0.0001:
        print("The total intensity is not 1, please verify the isotope list!")
        continue

    relAbandunceList = []
    for i in intensityList:
        relAbandunceList.append(100 * i / maxIntensity)

    # plot the spectrum

    print("Plotting spectrum...")

    plt.xlim(aveMass - massRange / 2, aveMass + massRange / 2)
    plt.ylim(0, 100)
    plt.xlabel('Mass (Da)')
    plt.ylabel('Relative Abundance')
    plt.title('Calculated isotope distribution of ' + compProp[0])
    plt.vlines(massList, [0], relAbandunceList, colors='r')
    plt.grid(True)

    print("Process completed!")

    plt.show()