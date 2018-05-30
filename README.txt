Below is a brief description of the bioinformatic pipeline used for
Pope et al. mBio 2017, Figure 8.

A new code repository was not created for this paper because it is not warranted.


1. Compute gene content dissimilarity for Actinobacteriophage_789 database:

    Script: pham_analysis.py in genome_analysis repo (tagged as v2.0)
    
    Input files:
        a. actino789 gene data (generated in Excel from SQL query)
        b. actino789 phage names (generated in Excel from SQL query)



2. Compute gene content dissimilarity for Cyanobacteriophage_209 database:

    Script: pham_analysis.py in genome_analysis repo (tagged as v2.0)
    
    Input files:
        a. cyano209 gene data (generated in Excel from SQL query)
        b. cyano209 phage names (generated in Excel from SQL query)



3. Compute maximum gene content dissimilarity gap for all phages in 
    Actinobacteriophage_789 and Cyanobacteriophage_209 databases:
    
    Script: analyze_mash_network.py in genome_analysis repo (tagged as v2.0)
    
    Input files:
        a. actino789 and cyano_209 phage names (generated in Excel from SQL query)
        b. actino789 and cyano_209 pham data (generated in pope_mbio_2017_analysis.R
            script)



4. Generate graphs for Figure 8:

    Script: pope_mbio_2017_analysis.R
    
    Input files:
        a. actino789 unique pairwise comparison list
            (generated in Excel)

        b. actino789 pairwise gene content dissimilarity
            (generated from pham_analysis.py)

        c. actino789 phage metadata
            (generated in Excel from SQL query)

        d. cyano209 unique pairwise comparison list
            (generated in Excel)

        e. cyano209 pairwise gene content dissimilarity
            (generated from pham_analysis.py)

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

        
        
