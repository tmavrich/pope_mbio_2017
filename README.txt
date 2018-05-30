Below is a brief description of the bioinformatic pipeline used for
Pope et al. mBio 2017, Figure 8.
The indicated scripts can be found in the pope_mbio_2017 repository.

Note: processing of the data using the scripts is not completely linear. Output from one script
(such as analyze_pham_data.py) gets imported into R where analyze_all_data.R processes the
data and outputs a new file for another script (such as analyze_mash_network.py), and
the new output gets imported back into R and processed by other steps in the
analyze_all_data.R code. Also, some input files are prepared manually using Excel with
data retrieved from the phage SQL database.



1. Compute gene content dissimilarity for Actinobacteriophage_789 database:

    Script: analyze_pham_data.py
    
    Input files:
        a. actino789 gene data (generated in Excel from SQL query)
        b. actino789 phage names (generated in Excel from SQL query)



2. Compute gene content dissimilarity for Cyanobacteriophage_209 database:

    Script: analyze_pham_data.py
    
    Input files:
        a. cyano209 gene data (generated in Excel from SQL query)
        b. cyano209 phage names (generated in Excel from SQL query)



3. Compute maximum gene content dissimilarity gap for all phages in 
    Actinobacteriophage_789 and Cyanobacteriophage_209 databases:
    
    Script: analyze_mash_network.py (note: this scripts receives gene content and 
    nucleotide distance data. If only using for gene content dissimilarity analysis, 
    nucleotide distance data can be set to 0.)
    
    Input files:
        a. actino789 and cyano_209 phage names (generated in Excel from SQL query)
        b. actino789 and cyano_209 pham data (generated in analyze_all_data.R
            script)



4. Generate graphs for Figure 8:

    Script: analyze_all_data.R
    
    Input files:
        a. actino789 unique pairwise comparison list
            (generated in Excel)

        b. actino789 pairwise gene content dissimilarity
            (generated from analyze_pham_data.py)

        c. actino789 phage metadata
            (generated in Excel from SQL query)

        d. cyano209 unique pairwise comparison list
            (generated in Excel)

        e. cyano209 pairwise gene content dissimilarity
            (generated from analyze_pham_data.py)

        f. cyano209 phage metadata
            (generated in Excel from SQL query)

        g. gcd gap analysis (all phages)
            (generated in Excel from analyze_mash_network.py output)

        h. gcd gap analysis (specific clusters)
            (generated in Excel from analyze_mash_network.py output)

        i. gcd gap analysis (specific subclusters)
            (generated in Excel from analyze_mash_network.py output)

        j. gcd gap analysis (specific host genera)
            (generated in Excel from analyze_mash_network.py output)

        
        
