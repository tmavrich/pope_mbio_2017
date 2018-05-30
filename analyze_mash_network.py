#Python3 script to perform misc analyses on mash networks
#Travis Mavrich
#20161206


import csv, os, sys, time

#Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:
    processed_file = sys.argv[1]
    genome_file = sys.argv[2]

except:

    print("\n\n\
        This is a Python3 script to perform misc analyses on phage networks\n\
        Script requires two arguments:\n\n\
        Execute script in any working directory\n\
        First argument: input file, which is the processed/filtered Mash data (csv-formatted)\n\
            0 = reference\n\
            1 = query\n\
            2 = nucleotide distance\n\
            3 = gene content dissimilarity\n\
            4+ = ...\n\
            \n\
        Second argument: input file of all phage names used in the Mash analysis (csv-formatted)\n\
            0 Phage Name (NOT phageID)\n\
            1 Predetermined group designation (Cluster, Subcluster, par phages, etc.)\n\n\
            2+ = ...\n\
        \n\n\n\
        The script outputs three file(s):\n\
        First output file: network analysis output (csv-formatted)\n\
            0 = group\n\
            1 = group size\n\
            2 = inter-group analysis: unique intra-group phage tally\n\
            3 = inter-group analysis: total # of intra-group phages with inter-group connections\n\
            4 = inter-group analysis: unique inter-group phage tally  (not including Unspecified groups)\n\
            5 = inter-group analysis: total # of inter-group phages that are connected to this particular group  (not including Unspecified groups)\n\
            6 = inter-group analysis: unique groups connected to this particular group tally (not including Unspecified groups)\n\
            7 = intra-group analysis: unique intra-group phage tally\n\
            8 = intra-group analysis: total # of intra-group phages with intra-group connections\n\
        Second output file: evolutionary mode analysis (csv-formatted)\n\
            0 = phage\n\
            1 = hgcf tally\n\
            2 = lgcf_tally\n\
            3 = out of range tally\n\
            4 = hgcf percent\n\
            5 = lgcf percent\n\
            6 = mode\n\
        Third output file: gene content dissimilarity gap analysis (csv-formatted)\n\
            0 = phage\n\
            1 = maximum gene content dissimilarity gap\n\
        \n")



    sys.exit(1)




#Expand home and working directory
home_dir = os.path.expanduser('~')
working_dir = os.path.abspath('.')




#Verify the genome data file exists

#Expand the path if it references the home directory
if genome_file[0] == "~":
    genome_file = home_dir + genome_file[1:]

if genome_file[0] != "/":
    genome_file = working_dir + '/' + genome_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
genome_file = os.path.abspath(genome_file)


if os.path.exists(genome_file) == False:
    print("\n\nInvalid genome data file path.\n\n")
    sys.exit(1)





#Verify the network data file exists

#Expand the path if it references the home directory
if processed_file[0] == "~":
    processed_file = home_dir + processed_file[1:]

if processed_file[0] != "/":
    processed_file = working_dir + '/' + processed_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
processed_file = os.path.abspath(processed_file)

if os.path.exists(processed_file) == False:
    print("\n\nInvalid network data file path.\n\n")
    sys.exit(1)





start_time = time.strftime("%x %X")





#Step 1: Import raw data
print("Importing data")

file_object = open(genome_file,"r")
file_reader = csv.reader(file_object)
genome_data_list = []
for row in file_reader:
    if row[1] == '\\N':
        row[1] = 'Singleton'
    genome_data_list.append(row)
file_object.close()



#Open the processed mash file
processed_data_list = []
processed_data_handle = open(processed_file, 'r')
processed_reader = csv.reader(processed_data_handle,delimiter=',')
for row in processed_reader:
    processed_data_list.append(row)
processed_data_handle.close





#Create output directory
date = time.strftime("%Y%m%d")
output_filename = '%s_network_mode_analysis' % date
output_folder = output_filename

try:
    os.mkdir(output_folder)
except:
    print("\nUnable to create output folder: %s" % output_folder)
    sys.exit(1)

os.chdir(output_folder)







#Create a set of phage groups
print("Creating sorted group list")
group_set = set()
for genome_data in genome_data_list:
    group_set.add(genome_data[1])
sorted_group_list = list(group_set)
sorted_group_list.sort()

#Create a group-phage dictionary.
#Key: group
#Value: set of phages
print("Creating group-phage dictionary")
group_phage_dict = {}
for group_id in sorted_group_list:
    matching_phage_set = set()
    for phage_group_data in genome_data_list:
        if phage_group_data[1] == group_id:
            matching_phage_set.add(phage_group_data[0])
    group_phage_dict[group_id] = matching_phage_set

#Create a phage-group dictionary.
#Key: phage
#Value: group
print("Creating phage-group dictionary")
phage_group_dict = {}
for genome_data in genome_data_list:
    phage_group_dict[genome_data[0]] = genome_data[1]

#Create a group-results dictionary.
#Key: group
#Value: network analysis results
#0 = inter-group analysis: unique intra-group phages
#1 = inter-group analysis: total # of intra-group phages with inter-group connections
#2 = inter-group analysis: unique inter-group phages
#3 = inter-group analysis: total # of inter-group phages that are connected to this particular group
#4 = inter-group analysis: unique groups connected to this particular group
#5 = intra-group analysis: unique intra-group phages
#6 = intra-group analysis: total # of intra-group phages with inter-group connections
print("Creating group-results dictionary")
group_results_dict = {}
for group_id in sorted_group_list:
    group_results_dict[group_id] = [set(),0,set(),0,set(),set(),0]



#Prepare results file
results_file_handle = open(str(output_filename + '.csv'),"w")
results_file_writer = csv.writer(results_file_handle)
column_headers = ["Group",\
                    "Group size",\
                    "Inter-group analysis: Intra-group unique phage tally",\
                    "Inter-group analysis:Intra-group total connections",\
                    "Inter-group analysis:Inter-group unique phage tally",\
                    "Inter-group analysis:Inter-group total connections",\
                    "Inter-group analysis:Inter-group unique group tally",\
                    "Intra-group analysis: Intra-group unique phage tally",\
                    "Intra-group analysis: Intra-group total connections"]
results_file_writer.writerow(column_headers)











print("Network analysis...\n\n")
for line in processed_data_list:


    ref_group = phage_group_dict[line[0]]
    query_group = phage_group_dict[line[1]]


    #If the group is "Unspecified", then skip. The analysis should only take into account specified Group data, with no missing values
    if ref_group == 'Unspecified' or query_group == 'Unspecified':
        continue


    ref_results = group_results_dict[ref_group]
    query_results = group_results_dict[query_group]


    #If the comparison is between phages of the same group, update intra-group analysis fields
    if ref_group == query_group:


        ref_results[5].add(line[0])
        ref_results[6] += 1

        query_results[5].add(line[1])
        query_results[6] += 1


    #If the comparison is between phages of different groups, update inter-group analysis fields
    else:


        ref_results[0].add(line[0])
        ref_results[1] += 1
        ref_results[2].add(line[1])
        ref_results[3] += 1
        ref_results[4].add(query_group)

        query_results[0].add(line[1])
        query_results[1] += 1
        query_results[2].add(line[0])
        query_results[3] += 1
        query_results[4].add(ref_group)






#Output all results to file
for key in group_results_dict:

    output_format = group_results_dict[key]
    output_list = [key,\
                    len(group_phage_dict[key]),\
                    len(output_format[0]),\
                    output_format[1],\
                    len(output_format[2]),\
                    output_format[3],\
                    len(output_format[4]),\
                    len(output_format[5]),\
                    output_format[6]]

    results_file_writer.writerow(output_list)



#Close the file
results_file_handle.close()








###
#Analyze mash and gene content dissimilarity patterns
print("Evolutionary mode analysis...\n\n")


#Create a set of all phages in the analysis data
analysis_phage_set = set()
for line in processed_data_list:
    analysis_phage_set.add(line[0])
    analysis_phage_set.add(line[1])


#Create a phage-data dictionary
#Key = phage
#Value = list of nucleotide distance and gene content dissimilarity data
print("Creating phage-data dictionary")
phage_data_dict = {}
phage_count = 0
for phage in analysis_phage_set:
    phage_count += 1
    print("Phage %s: %s" %(phage_count,phage))
    data_list = []
    for line in processed_data_list:
        if phage == line[0] or phage == line[1]:
            data_list.append([float(line[2]),float(line[3])])
    phage_data_dict[phage] = data_list




#Now compute the evolutionary mode for each phage
print("Computing evolutionary mode")
mode_data_list = []
for phage in phage_data_dict.keys():

    data_list = phage_data_dict[phage]
    hgcf_tally = 0
    lgcf_tally = 0
    out_of_range_tally = 0


    #Iterate through each set of coordinates
    for coordinates in data_list:


        #First determine if the coordinates fall within the measurable range to identify mode

        if (coordinates[0] < 0.06 and coordinates[1] < 0.22) or \
            (coordinates[0] > 0.28 and coordinates[1] > 0.79):
           out_of_range_tally += 1


        elif coordinates[0] < 0.16:

            if coordinates[1] > (3.5 * coordinates[0]):
                hgcf_tally += 1
            else:
                lgcf_tally += 1
        else:

            if coordinates[1] > ((2 * coordinates[0]) + 0.25):
                hgcf_tally += 1
            else:
                lgcf_tally += 1


    #Now compute the mode
    total_tally = hgcf_tally + lgcf_tally + out_of_range_tally
    in_range_tally = hgcf_tally + lgcf_tally

    #Make sure there is at least one measurable data point in either category
    if in_range_tally > 0:
        hgcf_percent = round(hgcf_tally/in_range_tally,2)
        lgcf_percent = round(lgcf_tally/in_range_tally,2)


        if hgcf_percent == 1:
            mode_exact = "hgcf"
        elif hgcf_percent == 0:
            mode_exact = "lgcf"
        else:
            mode_exact = "mixed"


        if hgcf_percent > 0.8:
            mode_80 = "hgcf"
        elif hgcf_percent < 0.2:
            mode_80 = "lgcf"
        else:
            mode_80 = "mixed"

    else:
        hgcf_percent = 0
        lgcf_percent = 0
        mode_exact = "unknown"
        mode_80 = "unknown"



    mode_data_list.append([phage,\
                            hgcf_tally,\
                            lgcf_tally,\
                            out_of_range_tally,\
                            hgcf_percent,\
                            lgcf_percent,\
                            mode_exact,\
                            mode_80])



#After all mode data is computed, output to file
mode_filename = '%s_mode_prediction' % date
mode_file_handle = open(str(mode_filename + '.csv'),"w")
mode_file_writer = csv.writer(mode_file_handle)
column_headers = ['phage',\
                'hgcf_tally',\
                'lgcf_tally',\
                'out_of_range_tally',\
                'hgcf_percent',\
                'lgcf_percent',\
                'mode_exact',\
                'mode_approx_80_percent']

mode_file_writer.writerow(column_headers)

for data in mode_data_list:
    mode_file_writer.writerow(data)



#Close the file
mode_file_handle.close()









#Now compute isolation or discreteness for each phage
#Use variables already computed for evolutionary mode
print("Computing phage discreteness")

gcd_gap_output_list = []

for phage in phage_data_dict:

    #print(phage)
    #input()
    data_list = phage_data_dict[phage]

    #print(data_list)
    #input()


    gcd_data_list = [0] #Add a zero fake data point, representing the phage against itself
    gcd_gap_list = []


    #Create list of only gene content dissimilarity data
    for element in data_list:
        #print(element)
        #input()
        gcd_data_list.append(element[1])

    #print(gcd_data_list)
    gcd_data_list.sort()
    #print(gcd_data_list)
    #input()
    index = 0
    while index < (len(gcd_data_list) - 1):
        gap = gcd_data_list[index+1] - gcd_data_list[index]
        gcd_gap_list.append(gap)
        index += 1

    gcd_gap_max = max(gcd_gap_list)
    gcd_gap_output_list.append([phage,gcd_gap_max])



#After all mode data is computed, output to file
gcd_gap_filename = '%s_gcd_gap_analysis' % date
gcd_gap_file_handle = open(str(gcd_gap_filename + '.csv'),"w")
gcd_gap_file_writer = csv.writer(gcd_gap_file_handle)
column_headers = ['phage',\
                'gcd_gap_max']

gcd_gap_file_writer.writerow(column_headers)
for element in gcd_gap_output_list:
    gcd_gap_file_writer.writerow(element)


#Close the file
gcd_gap_file_handle.close()










#Exit script
end_time = time.strftime("%x %X")

print("Start time: %s" %start_time)
print("Stop time: %s" %end_time)


print("Network and evolutionary mode analysis script complete.")
