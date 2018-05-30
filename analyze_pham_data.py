#Python3 script to compute conservation of phams within pre-determined groups of phages and/or
#proportion of shared phams between every pair of phages

#General notes:
#1. Script does not account for genes that straddle genome ends when calculating gene length. So it will not properly handle this edge case.
#2. For pham function analysis, all phams found in gene data are expected to be present in the pham function list, and vice versa. If this
#   is not the case, the script will exit.
#3. Pham proportion analysis outputs redundant data, so for two phages, it outputs (phage1, phage2, proportion data) as well as (phage2, phage1, proportion data)
#4. Gene data for phages that are not in the phage list does not cause errors. However, their associated pham data does get included in the analysis
#   (where their description gets added to that pham's dictionary). In general though, it is a good idea to ensure the phage, gene, and function input files
#   are consistent.

#Import requisite modules
import time, sys, os, csv, statistics


#Verify the correct arguments are provided, otherwise print description of the script and its required arguments.
try:
    genome_file = sys.argv[1]
    gene_file = sys.argv[2]
    function_file = sys.argv[3]

except:


    print("\n\n\
        This is a Python3 script to compute pham conservation frequencies and shared pham proportions\n\
        Execute script in any working directory\n\
        Script requires three input files:\n\n\
        First file: genome file (csv-formatted)\n\
            0 Phage Name\n\
            1 Predetermined group designation (Cluster, Subcluster, etc.)\n\n\
            \n\
        Second file: gene file (csv-formatted)\n\
            0 Phage Name\n\
            1 GeneID\n\
            2 Gene Start (1-based indexing)\n\
            3 Gene Stop (1-based indexing)\n\
            4 Pham\n\
            5 Gene description\n\
            6 Gene length\n\
            \n\
            \n\
        Third file: pham function file (csv-formatted)\n\
            0 Pham\n\
            1 Function\n\
        (Note: if no pham functions are available, enter 'none' for this argument)\n\n\n\n\n\n\
        Depending on the analyses performed, the following files are generated:\n\
        \n\
        1. Pham dictionary (Option 0) = list of all gene descriptions associated with each pham (csv-formatted)\n\
            0 pham\n\
            1 description\n\
        2. Genome-specific function summary (Option 0 with function analysis) = analysis of how functions contribute to gene totals and gene coding sequence\n\
            0 Phage Name\n\
            1 Total gene count\n\
            2...Function-specific gene count\n\
            3...Function-specific nucleotide gene length total\n\
            4...Function-specific nucleotide gene length mean\n\
            ....Repeat function-specific columns for each function\n\
            Last...Processed gene count\n\
        3. Group-specific pham tables (Option 0) = Pham presence/absence tables for each group\n\
            0 Pham\n\
            1+ All phages in group\n\
        4. Total pham conservation (Option 1) = list of all phams in the dataset and pham conservation at the Group level (csv-formatted)\n\
            0 Pham\n\
            1 Pham descriptions\n\
            2+...Each unique element/name in the Group designation\n\
        5. Genome-specific pham conservation (Option 1) = unique file for each phage, reflecting levels of pham conservation at the Group level (csv-formatted)\n\
            0 Phage Name\n\
            1 GeneID\n\
            2 Gene Start\n\
            3 Gene Stop\n\
            4 Pham\n\
            5 Gene description\n\
            6 Gene length\n\
            7 Pham descriptions\n\
            8+...Each unique element/name in the Group designation\n\
        6. Pairwise pham proportions (Option 2) = proportion of phams shared between each pair of phages (csv-formatted)\n\
            0 Phage1 name\n\
            1 Phage1 number of unshared phams\n\
            2 Phage1 shared proportion\n\
            3 Phage2 name\n\
            4 Phage2 number of unshared phams\n\
            5 Phage2 shared proportion\n\
            6 Number of shared phams\n\
            7 Average_shared_proportion\n\
            8 Jaccard similarity\n\
            9 Shared pham distribution mean\n\
            10 Shared pham distribution median\n\
            11 Shared pham distribution max\n\
            12 Unshared pham distribution mean\n\
            13 Unshared pham distribution median\n\
            14 Unshared pham distribution max\n\
            15 Unshared orpham count\n\
            16...Function-specific phage1 number of unshared phams\n\
            17...Function-specific phage2 number of unshared phams\n\
            18...Function-specific number of shared phams\n\
            19...Function-specific average shared proportion\n\
            ....Repeat function-specific columns for each function\n\
        7. Group connectedness (Option 3) = number of phams per group that are found in other groups or not (csv-formatted)\n\
            0 Group\n\
            1 Number of genomes\n\
            2 Total number of phams\n\
            3 Number of shared phams\n\
            4 Number of unshared phams\n\
            5 Mean number of groups in which the shared phams are present\n\
            6 Median number of groups in which the shared phams are present\n\
            7 Maximum number of groups in which the shared phams are present\n\
            8 Number of unshared phams present in all phages of the group\n\
            9+...Number of shared phams between each unique element/name in the Group designation\n\
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





#Verify the gene data file exists

#Expand the path if it references the home directory
if gene_file[0] == "~":
    gene_file = home_dir + gene_file[1:]

if gene_file[0] != "/":
    gene_file = working_dir + '/' + gene_file

#Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
gene_file = os.path.abspath(gene_file)

if os.path.exists(gene_file) == False:
    print("\n\nInvalid gene data file path.\n\n")
    sys.exit(1)



#Verify the pham function file exists

#Expand the path if it references the home directory
if function_file.lower() != "none":
    function_analysis = "yes"

    if function_file[0] == "~":
        function_file = home_dir + function_file[1:]

    if function_file[0] != "/":
        function_file = working_dir + '/' + function_file

    #Expand the path, to make sure it is a complete file path (in case user inputted path with './path/to/folder')
    function_file = os.path.abspath(function_file)

    if os.path.exists(function_file) == False:
        print("\n\nInvalid pham function file path.\n\n")
        sys.exit(1)
else:
    function_analysis = "no"








start_time = time.strftime("%x %X")
date = time.strftime("%Y%m%d")



print("\nWhich type of analysis do you want?\n\n")

print("0: Create pham dictionary, pham tables, and pham function summary.\n\
    No other pairwise or group comparisons are computed.\n\
    This is performed along with all other options.\n\n\n")

print("1: Pham conservation: For each unique pham for each individual genome,\n\
    the percentage of phages in each designated group (from the genome file)\n\
    that contain at least one copy of that pham is computed.\n\
    Therefore, a single output file is generated for each genome.\n\n\n")

print("2: Pairwise shared pham proportions: For each pairwise comparison of two genomes,\n\
    the proportion of shared phams between the two genomes\n\
    is computed and reported as a general or jaccard similarity.\n\n\n")

print("3: Shared phams between groups: For each unique group in the genome file,\n\
    the number of unique phams not found in any other group\n\
    and the number of phams found in at least one other group are computed.\n\n\n")

print("4: All\n")

valid_response = False
while valid_response == False:
    analysis_type = int(input("Choose analysis: "))
    if analysis_type >= 0 and analysis_type < 5:
        valid_response = True
    else:
        print("Invalid response.\n")



###Import raw data
print("Importing data")
database = input("Database name: ")

file_object = open(genome_file,"r")
file_reader = csv.reader(file_object)
genome_data_list = []
for row in file_reader:
    if row[1] == '\\N':
        row[1] = 'Singleton'
    genome_data_list.append(row)
file_object.close()


file_object = open(gene_file,"r")
file_reader = csv.reader(file_object)
gene_data_list = []
for row in file_reader:
    row[2] = int(row[2])
    row[3] = int(row[3])
    row[6] = int(row[6])
    gene_data_list.append(row)
file_object.close()




if function_analysis == "yes":
    file_object = open(function_file,"r")
    file_reader = csv.reader(file_object)
    function_data_list = []
    for row in file_reader:
        function_data_list.append(row)
    file_object.close()







#Create output directory to store all analyses
new_dir = date + "_" + database + "_pham_analysis"
os.mkdir(new_dir)
os.chdir(new_dir)






###Create dictionaries, sets, and lists to perform computations


#Create a set of phage genomes based on the genome data file
print("Creating phage set")
phage_set = set()
for genome_data in genome_data_list:
    phage_set.add(genome_data[0])
print("There are %s phage(s) in the phage set." %len(phage_set))



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


#Create pham set
print("Creating pham set")
pham_set = set()
for gene_data in gene_data_list:
    pham_set.add(gene_data[4])

#Create pham_description dictionary
#Key = pham
#Value = list of descriptions
print("Creating pham-description dictionary")
pham_description_dict = {}
for pham_id in pham_set:
    matching_description_list = []
    for gene_data in gene_data_list:
        if gene_data[4] == pham_id and gene_data[5] != '':
            matching_description_list.append(gene_data[5])
    matching_description_string = '; '.join(matching_description_list)
    pham_description_dict[pham_id] = matching_description_string




print("Exporting pham dictionary to csv file")
columns = ['pham','description']
file_name =  date + "_" + database + '_pham_dictionary.csv'
csvfile = open(file_name,'w')
writer = csv.writer(csvfile)
writer.writerow(columns)
for key in pham_description_dict:
    writer.writerow([key,pham_description_dict[key]])
csvfile.close()









#Create phage-gene dictionary
#Key = phage name
#Value = list of gene data
#0 = Phage name
#1 = GeneID
#2 = Gene start
#3 = Gene stop
#4 = Pham
#5 = Gene Description
#6 = Gene size
print("Creating phage-gene dictionary")
phage_gene_dict = {}
for genome_data in genome_data_list:
    matching_gene_list = []
    for gene_data in gene_data_list:
        if gene_data[0] == genome_data[0]:
            matching_gene_list.append(gene_data)
    phage_gene_dict[genome_data[0]] = matching_gene_list


#Create phage-pham dictionary
#Key = phage name
#Value = set of phams
print("Creating phage-pham dictionary")
phage_pham_dict = {}
for key in phage_gene_dict:
    matching_pham_set = set()
    for gene_data in phage_gene_dict[key]:
        matching_pham_set.add(gene_data[4])
    phage_pham_dict[key] = matching_pham_set




#Create group-pham dictionary
#Key = group
#Value = set of all phams found in the group
group_pham_dict = {}
print("Creating group-pham dictionary")
for group in group_phage_dict:
    group_pham_set = set()
    for phage in group_phage_dict[group]:
        group_pham_set = group_pham_set | phage_pham_dict[phage]
    group_pham_dict[group] = group_pham_set



#Create pham-group dictionary
#Key = pham
#Value = set of all groups that contain at least one gene of this pham
pham_group_dict = {}
print("Creating pham-group dictionary")
for pham in pham_set:
    pham_group_set = set()
    for group in group_pham_dict:
        if pham in group_pham_dict[group]:
            pham_group_set.add(group)
    pham_group_dict[pham] = pham_group_set





#Create pham-phage dictionary
#Key = pham
#Value = set of all phages from the genome data file that contain at least one gene in that pham
pham_phage_dict = {}
print("Creating pham-phage dictionary")
for pham in pham_group_dict.keys():
    pham_phage_set = set()
    for phage in phage_set:
        if pham in phage_pham_dict[phage]:
            pham_phage_set.add(phage)
    pham_phage_dict[pham] = pham_phage_set





if function_analysis == "yes":
    #Create pham-function dictionary
    #Key = pham
    #Value = function
    pham_func_dict = {}
    function_set = set()
    print("Creating pham-function dictionary")
    for row in function_data_list:
        pham_func_dict[row[0]] = row[1]
        function_set.add(row[1])


    sorted_function_list = list(function_set)
    sorted_function_list.sort()

    #Create function-pham dictionary
    #Key = function
    #Value = list of phams with that assigned function
    func_pham_dict = {}
    print("Creating function-pham dictionary")
    for function in function_set:
        associated_pham_set = set()
        for pham in pham_func_dict.keys():
            if function == pham_func_dict[pham]:
                associated_pham_set.add(pham)
        func_pham_dict[function] = associated_pham_set


    #Compare pham lists
    pham_dataset_diff = pham_description_dict.keys() ^ pham_func_dict.keys()
    if len(pham_dataset_diff) > 0:
        print("The phams in the gene data file do not exactly match the phams in the function file.")
        print("There are %s phams that are not found in both sets." % len(pham_dataset_diff))
        print("Unable to proceed. Exiting script.")
        sys.exit(1)


    #Compute the proportion of the genome sequence genes of each function contribute to
    print("Creating pham function summary")
    dataset_function_summary = []

    for phage in phage_set:

        #Retrieve all genes in this phage
        genes_to_analyze = phage_gene_dict[phage]

        phage_function_summary = [phage,len(genes_to_analyze)]

        processed_gene_tally = 0
        #Iterate through each function
        for function in sorted_function_list:

            function_nuc_length_list = []

            #Retrieve set of all phams with this function
            phams_in_function = func_pham_dict[function]


            for gene_data in genes_to_analyze:

                if gene_data[4] in phams_in_function:

                    #Compute length. Add +1 since 1-based indexing
                    function_nuc_length_list.append(abs(gene_data[3]-gene_data[2])+1)

            function_gene_count = len(function_nuc_length_list)
            function_nuc_gene_length_total = sum(function_nuc_length_list)

            if len(function_nuc_length_list) > 0:
                function_nuc_gene_length_mean = round(statistics.mean(function_nuc_length_list),2)
            else:
                function_nuc_gene_length_mean = 0

            phage_function_summary.extend([function_gene_count,function_nuc_gene_length_total,function_nuc_gene_length_mean])
            processed_gene_tally += function_gene_count

        phage_function_summary.append(processed_gene_tally)
        dataset_function_summary.append(phage_function_summary)

    print("Exporting function summary to csv file")

    header = []
    header.append(['Phage function summary'])
    header.append([new_dir])

    columns = ['phage',\
                'total_gene_count']

    #Add the function-specific column names
    for function in sorted_function_list:
        columns.append(function + "_gene_count")
        columns.append(function + "_nuc_gene_length_total")
        columns.append(function + "_nuc_gene_length_mean")

    columns.append('function_processed_gene_count')

    file_name =  new_dir + '_phage_function_summary.csv'
    csvfile = open(file_name,'w')
    writer = csv.writer(csvfile)
    for row in header:
        writer.writerow(row)
    writer.writerow(columns)
    for line in dataset_function_summary:
        writer.writerow(line)
    csvfile.close()






#Create pham tables based on group designation
print("Creating pham tables")


#Create output directory and move there
new_dir = date + "_" + database + "_pham_tables"

os.mkdir(new_dir)
os.chdir(new_dir)

#Iterate through each group
for group in group_set:

    group_output_list = []

    #Retrieve all phams in group
    phams_in_group = group_pham_dict[group]
    phams_in_group_sorted = list(phams_in_group)
    phams_in_group_sorted.sort()

    #Retrieve all phages in group
    phages_in_group = group_phage_dict[group]
    phages_in_group_sorted = list(phages_in_group)
    phages_in_group_sorted.sort()


    #Iterate through every pham in the group
    for pham in phams_in_group_sorted:

        pham_output_list = [pham]

        #For every phage in the group, determine if it has that pham
        for phage in phages_in_group_sorted:

            #Retrieve the list of phams present in this phage
            phams_in_phage = phage_pham_dict[phage]

            if pham in phams_in_phage:
                pham_output_list.append(1)
            else:
                pham_output_list.append(0)

        group_output_list.append(pham_output_list)


    pham_table_filename = "%s_pham_table.csv" %group
    csvfile = open(pham_table_filename,'w')
    writer = csv.writer(csvfile)

    columns = ["pham"]
    columns.extend(phages_in_group_sorted)
    writer.writerow(columns)
    for line in group_output_list:
        writer.writerow(line)
    csvfile.close()



#Move back to original working directory
os.chdir('..')






#The following code gets executed based on which type of analysis was chosen
#Pham conservation
if (analysis_type == 1) or (analysis_type == 4):



    ###Compute frequency of pham inclusion per group
    #Key = pham
    #Value = list of conservation values, each element corresponding to each group
    print("Computing pham conservation within groups")
    pham_conservation_list_dict = {}
    for pham in pham_set:
        pham_conservation_list = []
        index = 0
        while index < len(sorted_group_list):
            phages_to_inspect = group_phage_dict[sorted_group_list[index]]
            phage_sum = len(phages_to_inspect) #Total number of phages in this group
            pham_sum = 0
            for phage in phages_to_inspect:
                if pham in phage_pham_dict[phage]:
                    pham_sum += 1
            pham_conservation_list.append(pham_sum/phage_sum)
            index += 1
        pham_conservation_list_dict[pham] = pham_conservation_list

    #Now output the entire conservation dataset
    total_conservation_header = []
    total_conservation_header.append(['Pham conservation script'])
    total_conservation_header.append(['All pham conservation data'])
    total_conservation_header.append(['Database: %s' %database])
    total_conservation_header.append([date])

    total_conservation_columns = ['Pham','PhamDescriptions']
    index4 = 0
    while index4 < len(sorted_group_list):
        total_conservation_columns.append(sorted_group_list[index4])
        index4 += 1

    file_name = '%s_%s_total_pham_conservation.csv' %(date,database)
    csvfile = open(file_name,'w',newline='')
    writer = csv.writer(csvfile)
    for row in total_conservation_header:
        writer.writerow(row)
    writer.writerow(total_conservation_columns)

    for pham in pham_conservation_list_dict:
        result_list = pham_conservation_list_dict[pham]
        pham_descriptions = pham_description_dict[pham]
        output_list = []
        output_list.append(pham)
        output_list.append(pham_descriptions)
        output_list.extend(result_list)
        writer.writerow(output_list)
    csvfile.close()









    ###Output conservation data for each phage
    #Creates a separate csv file for each phage. Each row is a gene, each column is a pham conservation per group


    #Create output directory and move there
    new_dir = date + "_" + database + "_pham_conservation"
    os.mkdir(new_dir)
    os.chdir(new_dir)


    print("Exporting data to csv files")
    header = []
    header.append(['Pham conservation script'])
    header.append([new_dir])

    columns = ['Phage','GeneID','GeneStart','GeneStop','Pham','GeneDescription','GeneSize','PhamDescriptions']
    index2 = 0
    while index2 < len(sorted_group_list):
        columns.append(sorted_group_list[index2])
        index2 += 1


    for group in sorted_group_list:
        for phage in group_phage_dict[group]:
            file_name = group + '_' + phage + '_pham_conservation.csv'
            csvfile = open(file_name,'w',newline='')
            writer = csv.writer(csvfile)
            for row in header:
                writer.writerow(row)
            file_name_list = [file_name]
            writer.writerow(file_name_list)
            writer.writerow(columns)

            for gene in phage_gene_dict[phage]:
                pham_descriptions = pham_description_dict[gene[4]]

                output_list = []
                output_list.append(gene[0])
                output_list.append(gene[1])
                output_list.append(gene[2])
                output_list.append(gene[3])
                output_list.append(gene[4])
                output_list.append(gene[5])
                output_list.append(gene[6])
                output_list.append(pham_descriptions)

                result_list = pham_conservation_list_dict[gene[4]]
                index3 = 0
                while index3 < len(result_list):
                    output_list.append(result_list[index3])
                    index3 += 1
                writer.writerow(output_list)
            csvfile.close()

    #Move back to original working directory
    os.chdir('..')



#Pairwise shared pham proportions
if (analysis_type == 2) or (analysis_type == 4):



    #Create output directory and move there
    new_dir = date + "_" + database + "_pairwise_pham_proportions"
    os.mkdir(new_dir)
    os.chdir(new_dir)

    ###Compute proportions of shared phams between every pair of phages
    print("Computing pairwise shared pham proportions")

    #This list will store the proportions of phams shared between phages
    proportion_data_list = []


    for phage1 in phage_pham_dict:

        phage1_pham_set = phage_pham_dict[phage1]

        #If there are no phams in the first genome, then skip this genome
        if len(phage1_pham_set) == 0:
            print("Phage %s contains no phams." % phage1)
            continue

        for phage2 in phage_pham_dict:

            #Reset these variables for every pairwise comparison
            proportion_data_summary = []
            shared_pham_distribution_list = [] # number of groups that contain >1 gene in this pham
            unshared_pham_distribution_list = [] # number of groups that contain >1 gene in this pham
            unshared_orpham_list = [] # number of phams in this comparison that are only present once in the entire dataset


            #If both phages to compare are the same, then skip this genome
            if phage1 == phage2:
                print("Phages %s and %s are the same." % (phage1,phage2))
                continue


            phage2_pham_set = phage_pham_dict[phage2]

            #If there are no phams in the second genome, then skip this genome
            if len(phage2_pham_set) == 0:
                print("Phage %s contains no phams." % phage2)
                continue


            intersection_set = phage1_pham_set & phage2_pham_set
            union_set = phage1_pham_set | phage2_pham_set


            phage1_unshared_pham_set = phage1_pham_set - phage2_pham_set
            phage2_unshared_pham_set = phage2_pham_set - phage1_pham_set

            phage1_unshared_size = len(phage1_unshared_pham_set)
            phage2_unshared_size = len(phage2_unshared_pham_set)
            phage1_shared_proportion = len(intersection_set)/len(phage1_pham_set)
            phage2_shared_proportion = len(intersection_set)/len(phage2_pham_set)
            ave_proportion = (phage1_shared_proportion + phage2_shared_proportion)/2
            jaccard_similarity = len(intersection_set)/len(union_set)




            #Now compute the group distribution of each shared and unshared pham
            for pham in intersection_set:
                shared_pham_distribution_list.append(len(pham_group_dict[pham]))

            for pham in phage1_unshared_pham_set:
                unshared_pham_distribution_list.append(len(pham_group_dict[pham]))

                if len(pham_phage_dict[pham]) == 1:
                    unshared_orpham_list.append(pham)

            for pham in phage2_unshared_pham_set:
                unshared_pham_distribution_list.append(len(pham_group_dict[pham]))

                if len(pham_phage_dict[pham]) == 1:
                    unshared_orpham_list.append(pham)


            if len(shared_pham_distribution_list) > 0:
                shared_pham_distribution_mean = statistics.mean(shared_pham_distribution_list)
                shared_pham_distribution_median = statistics.median(shared_pham_distribution_list)
                shared_pham_distribution_max = max(shared_pham_distribution_list)
            else:
                shared_pham_distribution_mean = 0
                shared_pham_distribution_median = 0
                shared_pham_distribution_max = 0


            if len(unshared_pham_distribution_list) > 0:
                unshared_pham_distribution_mean = statistics.mean(unshared_pham_distribution_list)
                unshared_pham_distribution_median = statistics.median(unshared_pham_distribution_list)
                unshared_pham_distribution_max = max(unshared_pham_distribution_list)
            else:
                unshared_pham_distribution_mean = 0
                unshared_pham_distribution_median = 0
                unshared_pham_distribution_max = 0


            proportion_data_summary.extend([phage1,\
                                        phage1_unshared_size,\
                                        round(phage1_shared_proportion,4),\
                                        phage2,\
                                        phage2_unshared_size,\
                                        round(phage2_shared_proportion,4),\
                                        len(intersection_set),\
                                        round(ave_proportion,4),\
                                        round(jaccard_similarity,4),\
                                        round(shared_pham_distribution_mean,4),\
                                        shared_pham_distribution_median,\
                                        shared_pham_distribution_max,\
                                        round(unshared_pham_distribution_mean,4),\
                                        unshared_pham_distribution_median,\
                                        unshared_pham_distribution_max,\
                                        len(unshared_orpham_list),\
                                        ])


            #Now compute pham function analysis, if selected
            if function_analysis == "yes":

                for function in sorted_function_list:


                    function_pham_set = func_pham_dict[function]
                    shared_function_pham_set = function_pham_set & intersection_set
                    phage1_unshared_function_pham_set = function_pham_set & phage1_unshared_pham_set
                    phage2_unshared_function_pham_set = function_pham_set & phage2_unshared_pham_set




                    #compute pham function-specific gene content dissimilarity
                    phage1_unshared_function_size = len(phage1_unshared_function_pham_set)
                    phage2_unshared_function_size = len(phage2_unshared_function_pham_set)
                    phage1_total_function_size = len(shared_function_pham_set) + phage1_unshared_function_size
                    phage2_total_function_size = len(shared_function_pham_set) + phage2_unshared_function_size


                    if phage1_total_function_size > 0:
                        phage1_shared_function_proportion = len(shared_function_pham_set)/phage1_total_function_size
                    else:
                        phage1_shared_function_proportion = -1


                    if phage2_total_function_size > 0:
                        phage2_shared_function_proportion = len(shared_function_pham_set)/phage2_total_function_size
                    else:
                        phage2_shared_function_proportion = -1


                    if phage1_shared_function_proportion == -1 and phage2_shared_function_proportion == -1:
                            ave_shared_function_proportion = -1

                    elif phage1_shared_function_proportion == -1 and phage2_shared_function_proportion != -1:
                            ave_shared_function_proportion = phage2_shared_function_proportion

                    elif phage1_shared_function_proportion != -1 and phage2_shared_function_proportion == -1:
                            ave_shared_function_proportion = phage1_shared_function_proportion

                    else:
                        ave_shared_function_proportion = (phage1_shared_function_proportion + phage2_shared_function_proportion)/2

                    proportion_data_summary.append(phage1_unshared_function_size)
                    proportion_data_summary.append(phage2_unshared_function_size)
                    proportion_data_summary.append(len(shared_function_pham_set))
                    proportion_data_summary.append(round(ave_shared_function_proportion,4))


            #After all data is computed, append to the master list of proportion data
            proportion_data_list.append(proportion_data_summary)


    ###Output shared pham proportion data and function analysis (if selected)
    print("Exporting proportion data to csv file")


    #Don't add a header to this file.
    #This file is large, and it gets loaded into R for downstream analysis, so the header rows complicates file import
    columns = ['phage1_name',\
                'phage1_number_of_unshared_phams',\
                'phage1_shared_proportion',\
                'phage2_name',\
                'phage2_number_of_unshared_phams',\
                'phage2_shared_proportion',\
                'number_of_shared_phams',\
                'average_shared_proportion',\
                'jaccard_similarity',\
                'shared_pham_distribution_mean',\
                'shared_pham_distribution_median',\
                'shared_pham_distribution_max',\
                'unshared_pham_distribution_mean',\
                'unshared_pham_distribution_median',\
                'unshared_pham_distribution_max',\
                'unshared_orpham_count',\
                ]


    #Add the function-specific column names if selected by user
    if function_analysis == "yes":
        print("Exporting function analysis to csv file")
        for function in sorted_function_list:
            columns.append(function + "_phage1_number_of_unshared_phams")
            columns.append(function + "_phage2_number_of_unshared_phams")
            columns.append(function + "_number_of_shared_phams")
            columns.append(function + "_average_shared_proportion")


    file_name =  new_dir + '.csv'
    csvfile = open(file_name,'w')
    writer = csv.writer(csvfile)
    writer.writerow(columns)
    for line in proportion_data_list:
        writer.writerow(line)
    csvfile.close()


    os.chdir('..')






#Compute group connectedness
if (analysis_type == 3) or (analysis_type == 4):

    #Make sure there is more than one group.
    if not len(sorted_group_list) > 1:
        print("Only one group in data set.")
        print("Unable to compute group connectedness.")
    else:
        #Create output directory and move there
        new_dir = date + "_" + database + "_group_pham_overlaps"
        os.mkdir(new_dir)
        os.chdir(new_dir)

        ###Compute proportions of phams within each group that are present in other groups
        print("Computing group connectedness")


        #This dictionary will store the set of phams in each group that are found in other groups
        #Key = group
        #Value = set of all phams of a group that are found at least once in other groups
        group_connections_dict = {}



        index5 = 0
        while index5 < len(sorted_group_list):


            group1 = sorted_group_list[index5]
            group1_total_pham_set = group_pham_dict[group1]
            group1_shared_pham_set = set() #How many phams in this group are found in at least one other group?
            group1_shared_pham_per_group_list = [] #Each element represent the number of shared phams between each particular group

            group1_genome_count = len(group_phage_dict[group1]) #How many phages in this group?
            group1_total_pham_count = len(group1_total_pham_set) #How many phams are in this group?

            index6 = 0
            while index6 < len(sorted_group_list):

                group2 = sorted_group_list[index6]
                #Compare pham sets of the two groups and compute the phams shared
                #Add the shared phams to the first group's main set of all phams shared in other groups
                if group1 != group2:

                    group_pham_intersection_set = group1_total_pham_set & group_pham_dict[group2]
                    group1_shared_pham_per_group_list.append(len(group_pham_intersection_set))
                    group1_shared_pham_set = group1_shared_pham_set | group_pham_intersection_set


                else:
                    group1_shared_pham_per_group_list.append(len(group1_total_pham_set))
                    print("Skipping same group comparison")

                index6 += 1


            #For the unshared phams, compute whether they are completely conserved in that group or not
            group1_unshared_pham_set = group1_total_pham_set - group1_shared_pham_set #How many phams in this group are not found in any other group?

            pham_group1_completely_conserved = 0
            for pham in group1_unshared_pham_set:

                pham_group_tally = 0
                for phage in group_phage_dict[group1]:

                    if pham in phage_pham_dict[phage]:
                        pham_group_tally += 1

                if pham_group_tally == len(group_phage_dict[group1]):
                    pham_group1_completely_conserved += 1


            #Now that the shared and unshared phams are computed, for each group compute the mean, median, and max number of groups each shared pham is found in
            group1_shared_pham_distribution = []
            for pham in group1_shared_pham_set:
                group1_shared_pham_distribution.append(len(pham_group_dict[pham]))

            if len(group1_shared_pham_distribution) > 0:
                group1_shared_pham_distribution_mean = statistics.mean(group1_shared_pham_distribution)
                group1_shared_pham_distribution_median = statistics.median(group1_shared_pham_distribution)
                group1_shared_pham_distribution_max = max(group1_shared_pham_distribution)
            else:
                group1_shared_pham_distribution_mean = 0
                group1_shared_pham_distribution_median = 0
                group1_shared_pham_distribution_max = 0


            #Now compile all data for each group in a list for output
            group1_summary_list = [group1,\
                                    group1_genome_count,\
                                    group1_total_pham_count,\
                                    len(group1_shared_pham_set),\
                                    len(group1_unshared_pham_set),\
                                    group1_shared_pham_distribution_mean,\
                                    group1_shared_pham_distribution_median,\
                                    group1_shared_pham_distribution_max,\
                                    pham_group1_completely_conserved]
            group1_summary_list.extend(group1_shared_pham_per_group_list)
            group_connections_dict[group1] = group1_summary_list



            #Move to the next index
            index5 += 1


        ###Output shared pham data
        print("Exporting group shared pham data to csv file")
        header = []
        header.append(['Number of phams per group that are found in other groups'])
        header.append([new_dir])

        columns = ['group','number_of_genomes',\
                'total_number_of_phams',\
                'number_of_shared_phams',\
                'number_of_unshared_phams',\
                'shared_pham_group_distribution_mean',\
                'shared_pham_group_distribution_median',\
                'shared_pham_group_distribution_max',\
                'number_of_unshared_phams_present_in_every_phage_of_the_group']
        #Add headers for each specific group in the group list
        columns.extend(sorted_group_list)


        file_name =  new_dir + '.csv'
        csvfile = open(file_name,'w')
        writer = csv.writer(csvfile)
        writer.writerow(columns)
        for group in group_connections_dict:
            writer.writerow(group_connections_dict[group])
        csvfile.close()

        os.chdir('..')








end_time = time.strftime("%x %X")

print("Start time: %s" %start_time)
print("Stop time: %s" %end_time)

print("\n\nPham analysis completed.")
