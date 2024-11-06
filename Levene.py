from scipy.stats import levene
from scipy.stats import ttest_ind
import pandas as pd

# load text files into an array

file_suffix = "_R1_density.txt"

chromosomes = [ "Tb427_01_v4", "Tb427_02_v4", "Tb427_03_v4", "Tb427_04_v4",
                "Tb427_05_v4", "Tb427_06_v4", "Tb427_07_v4", "Tb427_08_v4",
                "Tb427_09_v4", "Tb427_10_v5", "Tb427_11_01_v4"]

density_arrays = []
binary_arrays = []

for chromosome in chromosomes:
    file_name = f"Density_and_binary files/R1/Density/{chromosome}{file_suffix}"
    array = []
    with open(file_name, 'r') as file:
        for line in file:
            array.append(float(line.strip()))


    lower_bound = len(array)//10
    upper_bound = (len(array)//10) * 9
    array = array[lower_bound:upper_bound]
    
    density_arrays.append(array)
    file_name = f"Density_and_binary files/R1/Binary/{chromosome}_R1_is_ssr.txt"
    array = []
    with open(file_name, 'r') as file:
        for line in file:
            array.append(float(line.strip()))

    array = array[lower_bound:upper_bound]
    binary_arrays.append(array)

# Perform Levene's test
statistic, p_value = levene(*density_arrays)

print(f"Levene's test statistic: {statistic}")
print(f"p-value: {p_value}")
    
crosslink_density_tss = []
crosslink_density_non_tss = []

stats =[]

for i in range(len(chromosomes)):
    tss = []
    non_tss = []
    for j in range(len(density_arrays[i])):
        if binary_arrays[i][j] == 1:
            crosslink_density_tss.append(density_arrays[i][j])
            tss.append(density_arrays[i][j])
        else:
            crosslink_density_non_tss.append(density_arrays[i][j])
            non_tss.append(density_arrays[i][j])
    statistic, p_value = ttest_ind(tss, non_tss, equal_var=False)

    # print(f"t-test statistic for {chromosomes[i]}: {statistic}")
    print(f"p-value for {chromosomes[i]}: {p_value}")
    stats.append((chromosomes[i], statistic, p_value))


statistic, p_value = ttest_ind(crosslink_density_tss, crosslink_density_non_tss, equal_var=False)

print(f"t-test statistic for all chromosomes: {statistic}")
print(f"p-value for all chromosomes: {p_value}")

stats.append(("All chromosomes", statistic, p_value))

df = pd.DataFrame(stats, columns=["Chromosome", "Statistic", "P-value"])

df.to_csv("ttest_results_0.1_telomeric_cutoff.tsv", sep='\t', index=False)
df.to_excel("Ttest_results_0.1_telomeric_cutoff.xlsx", index=False)
