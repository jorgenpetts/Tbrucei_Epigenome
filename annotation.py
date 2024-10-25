import numpy as np

## TODO add code to write how far peaks are from the K125ac peak. and to identify different SSRs
# Define the file paths
file1 = 'RPB9_chIP.bam_peaks.narrowPeak'
file2 = 'H2AZ_chIP-seq.bam_peaks.narrowPeak'
file3 = 'BF1_K125ac_peaks.narrowPeak'

# Initialize empty arrays
RPB9 = np.array([])
H2AZ = np.array([])
K125ac = np.array([])

# Read data from file1

with open(file1, 'r') as f:
    for line in f:
        line = line.split('\t')
        chromName = line[0]
        chromStart = int(line[1])
        chromEnd = int(line[2])
        peakName = line[3]
        score = int(line[4])
        strand = line[5]
        signalValue = float(line[6])
        pValue = float(line[7])
        qValue = float(line[8])
        peaks = int(line[9]) + chromStart
        RPB9 = np.append(RPB9, [chromName, chromStart, chromEnd, peakName, score, strand, signalValue, pValue, qValue, peaks])
        # Process each line and append to array1
        # ...

RPB9 = RPB9.reshape(-1, 10)

# Read data from file2
with open(file2, 'r') as f:
    for line in f:
        line = line.split('\t')
        chromName = line[0]
        chromStart = int(line[1])
        chromEnd = int(line[2])
        peakName = line[3]
        score = int(line[4])
        strand = line[5]
        signalValue = float(line[6])
        pValue = float(line[7])
        qValue = float(line[8])
        peaks = int(line[9]) + chromStart
        H2AZ = np.append(H2AZ, [chromName, chromStart, chromEnd, peakName, score, strand, signalValue, pValue, qValue, peaks])
        # Process each line and append to array1
        # ...

H2AZ = H2AZ.reshape(-1, 10)
        

# Read data from file3
with open(file3, 'r') as f:
    for line in f:
        line = line.split('\t')
        chromName = line[0]
        chromStart = int(line[1])
        chromEnd = int(line[2])
        peakName = line[3]
        score = int(line[4])
        strand = line[5]
        signalValue = float(line[6])
        pValue = float(line[7])
        qValue = float(line[8])
        peaks = int(line[9]) + chromStart
        K125ac = np.append(K125ac, [chromName, chromStart, chromEnd, peakName, score, strand, signalValue, pValue, qValue, peaks])
        

K125ac = K125ac.reshape(-1, 10)



# Use 8kb upstream and downstream of the peak summit for K125ac
H2AZ_count = 0
RPB9_count = 0
average = 0
count = 0

tempH2AZ = np.array([])
tempRPB9 = np.array([])
final = np.array([])
final_annotation = np.array([])
for i in range(K125ac.shape[0]):
    name = K125ac[i, 0]
    start = K125ac[i, 1]
    end = K125ac[i, 2]
    peak = int(K125ac[i, 9])
    for j in range(H2AZ.shape[0]):
        if H2AZ[j, 0] == name:
            if peak - 8000 < int(H2AZ[j, 9]) < peak + 8000: # Check if H2AZ peak is within 8kb of K125ac peak
                H2AZ_count += 1
                tempH2AZ = np.append(tempH2AZ, H2AZ[j])

                #TODO Get all H2AZ peaks within 8kb of K125ac peak and store in temp, then do same for RPB9
                #TODO process the peaks and put into a new array that will be used for the next step to put into 
                #TODO .narowPeak file format
            if peak + 8000 < int(H2AZ[j, 9]):
                break
    
    for j in range(RPB9.shape[0]):
        if RPB9[j, 0] == name:
            if peak - 8000 < int(RPB9[j, 9]) < peak + 8000:
                RPB9_count += 1
                tempRPB9 = np.append(tempRPB9, RPB9[j])
            
            if peak + 8000 < int(RPB9[j, 9]):
                break
       

    if H2AZ_count != 0 and RPB9_count != 0:

        tempRPB9 = tempRPB9.reshape(-1, 10)
        tempH2AZ = tempH2AZ.reshape(-1, 10)
        if count == 0:
            print(f"K125ac Peak: {peak}")
        average += peak
        for j in range(tempH2AZ.shape[0]):
            if count == 0:
                print(f"Start: {tempH2AZ[j, 1]} Peak: {(tempH2AZ[j, 9])}")
            # add stuff for annotation here to final_annotation
            average += int(tempH2AZ[j, 9])
    
        for j in range(tempRPB9.shape[0]):
            if count == 0:
                print(f"Start: {tempRPB9[j, 1]} Peak: {(tempRPB9[j, 9])}")
            # add stuff for annotation here to final_annotation
            average += int(tempRPB9[j, 9])

        average = average // (tempH2AZ.shape[0] + tempRPB9.shape[0] + 1)
        if count == 0:
            print(average)
        count += 1
        finalName = f"K125ac_match_{count}"
        
        final = np.append(final, [name, average - 1000, average + 1000, finalName, K125ac[i, 4], '.', K125ac[i, 6], K125ac[i, 7], K125ac[i, 8], 1000])

    H2AZ_count = 0
    RPB9_count = 0
    average = 0
    tempH2AZ = np.array([])
    tempRPB9 = np.array([])

final = final.reshape(-1, 10)
print(final[0])

with open('annotatedPeaks.narrowPeak', 'w') as f:
    for row in final:
        f.write('\t'.join(map(str, row)) + '\n')

