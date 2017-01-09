#include <fstream>
#include <iostream>
#include <string>
#include <vector>
#include <sstream>
#include <iterator>
#include <cstdlib>
#include <algorithm>  
#include <cstring>
#include <map>
#include <cerrno>
#include <stdexcept>
#include <boost/iostreams/filtering_streambuf.hpp>
#include <boost/iostreams/copy.hpp>
#include <boost/iostreams/filter/gzip.hpp>
#include <zlib.h>


/*
To do : 
Add optional support for multple individuals
*/ 

/*
This program reads a vcf, a gvcf and plink map file and outputs a ped file subset 
down to the loci in the map file. This is useful to various programs in 
population genetics like Plink, Eigenstrat, King etc. 
which accept ped and map files as their input.

First the vcf is scanned and if the site is not found, it is looked up in the gvcf
file. For sites in the gvcf file if the coverage is < 10X,
a missing genotype is output. All sites found in the gvcf are output as 
homozygous reference.


gvcf : https://software.broadinstitute.org/gatk/guide/article?id=4017
map file : https://www.cog-genomics.org/plink2/formats#map
*/


void with_io_exceptions(std::ios& io)
{
    io.exceptions(std::ios_base::badbit | std::ios_base::failbit);
}

std::vector<std::string> split_line(const std::string &line, char delim)
/* Function to read a gvcf/vcf file , returns a vector of strings */
{
  // Use a boost solution for future  
    std::stringstream ss(line);
    std::string column;
    std::vector<std::string> elements;

    while (std::getline(ss,column,delim))
    {
	elements.push_back(column);
    }
    return elements;
}

std::string return_genotypes_on_dp(int dp, const std::string &ref, const std::string &alt)
/* Return genotypes based upon DP value
   This is for sites from the gvcf, if DP >= 10 
   store as homozygous reference
*/
{
    std::string temp_genotype;   
    if (dp >= 10)
    {
	temp_genotype = ref;
	temp_genotype += ' ';
	temp_genotype += ref;
    }
    else
	temp_genotype = "0 0"; // Set as missing genotype
    
    return temp_genotype;
}

std::string return_genotypes_on_gt(const std::string &gt, const std::string &ref, const std::string &alt)
/* Return genotypes based upon the GT
   field (zygocity)
*/
{
    std::string temp_genotypes;
    // GT is expected to be 3 characters , i.e. 1/1, 2/2, 1/2 etc.
    if (gt.length() != 3)
	throw "The GT field for the multiallelic site has not been resolved properly";
   
    if (gt == "0/0") // Homozygous ref
	temp_genotypes = ref + " " + ref;
    else if(gt == "1/1") // Homozygous alt
	temp_genotypes = alt + " " + alt;
    else if(gt == "1/0" or gt == "0/1" or gt == "1/2" or gt == "2/1") // Heterozygous site
	temp_genotypes = ref + " " + alt;
    else
	throw "The GT field for the multiallelic site has not been resolved properly";
    
    return temp_genotypes;
}

std::tuple< std::vector < std::vector<int> > , std::vector < std::vector<std::string> > > file_reader(const std::string& input_file, const std::string& file_type)
/* Read the gvcf and return the loci and genotypes as 2d arrays in a tuple */
{  
    //std::ifstream infile(input_file.c_str());
    std::ifstream file(input_file.c_str(),std::ios_base::in|std::ios_base::binary);
    boost::iostreams::filtering_streambuf<boost::iostreams::input> inbuf;
    inbuf.push(boost::iostreams::gzip_decompressor());
    inbuf.push(file);
    // Convert streambuf to istream
    std::istream infile(&inbuf);
         
    std::string line;    
    unsigned long int line_counter = 0;
    unsigned long int i = 0;
    int j;
    std::vector < std::vector<int> > positions(22,std::vector<int>()); // Store [i][j] -> positions at different loci 
    std::vector < std::vector<std::string> > genotype_mat(22,std::vector<std::string>()); // Store [i][j] -> genotypes at different loci
    std::string temp_genotypes;   
    std::string chrom;
    std::string prev_chrom;
    unsigned long int pos;
    std::string ref;
    std::string alt;
    std::string tmp_alt;
    unsigned long int end;
    int dp = 1; // Declare DP to 1 , for vcf we do not care about this value, will be overridden for gvcf
    std::string gt; // This will be used to decipher the zygocity of the site (1/1,1/0,etc.) , only used for vcf
    std::string tmp_gt;
    std::vector<std::string> temp_store;
    int info_len; // Length of the info field     
    std::vector<std::string> elements; // For storing the split contents of a gvcf/vcf line 
    std::vector<std::string> chroms{"1","2","3","4","5","6","7","8","9","10",
	    "11","12","13","14","15","16","17","18","19","20","21","22"};
  
    if (infile)
    {
	std::cout << "Reading the " << file_type << " file!\n";
	while (std::getline(infile,line))
	{
	    if (line[0] == '#')
		continue;
	    else
	    {
		elements = split_line(line,'\t');
		chrom = elements[0];
		if (std::find(chroms.begin(),chroms.end(),chrom) != chroms.end() == false)
		    continue;
		if (chrom != prev_chrom)
		    i = 0;
		
		pos = std::atoi(elements[1].c_str());
		if (elements[3].length() > 1) //
		    ref = "0";
		else
		    ref = elements[3][0];
		end = pos; // will overide this later if the file is a gvcf 

		// Need to process a bit further for getting alt allele
		temp_store = split_line(elements[4],',');
		
		// Get gvcf specific values i.e. END and DP
		if (file_type == "gvcf")
		{
		    if (temp_store.size() == 1)
		    {
			alt = temp_store[0][0];
			temp_store = split_line(elements[7],'=');
			end = std::atoi(temp_store[1].c_str());
		    }
		    			
		    //std::cout << end << " " << pos << '\n';
		    // Split the last field of the gvcf
		    temp_store = split_line(elements[9],':');
		    // Calculate the length of the split vector, dp position varies depending on the site being a variant or a ref site
		    info_len = sizeof(temp_store)/sizeof(temp_store[0]);
		    if (info_len == 6)
			dp = std::atoi(temp_store[2].c_str());
		    else
			dp = std::atoi(temp_store[1].c_str());

		    if (ref == "0")
			alt = "0";
		    temp_genotypes = return_genotypes_on_dp(dp,ref,alt);
		    int k2;
		    k2 = std::atoi(chrom.c_str()) - 1;
		    genotype_mat[k2].push_back(temp_genotypes);
		    positions[k2].push_back(pos);
		    ++line_counter;
		}
		
		else if(file_type == "vcf")
		{		    
		    tmp_alt = temp_store[0];
		    if (tmp_alt.length() > 1)
		    {
			ref = "0";
			alt = "0";
		    }
		    else
			alt = temp_store[0]; // Use only the first allele
		    
		    temp_store = split_line(elements[9],':');
		    info_len = sizeof(temp_store)/sizeof(temp_store[0]);
		    tmp_gt = temp_store[0];
		    gt = split_line(tmp_gt,',')[0]; // Use only the first allele
		    
		    dp = std::atoi(temp_store[2].c_str());		    
		    if (ref == "0")
			alt = "0";
		    
		    temp_genotypes = return_genotypes_on_gt(gt,ref,alt);
		    
		    int k2;
		    k2 = std::atoi(chrom.c_str()) - 1;
		    genotype_mat[k2].push_back(temp_genotypes);
		    positions[k2].push_back(pos);		    
		    ++line_counter;		    
		}		
		else
		    throw "Incorrect file type specified , specify either a gvcf or vcf file";
	    }
	}
    }
    else
       throw "Error Reading gvcf/vcf File!";
    
    //infile.close();
    file.close();
    std::cout << "Read : " << line_counter << " Lines!" << '\n';    
    return std::make_tuple(positions,genotype_mat);
}


std::vector< std::vector<int> > return_map_coordinates(const std::string& input_file)
{
    std::string line;
    std::ifstream infile(input_file.c_str());
    std::vector<std::string> elements;
    int chrom;
    unsigned long int pos;
    std::string temp;
    std::vector < std::vector<int> > positions(23,std::vector<int>());
    if(infile)
    {
	while (std::getline(infile,line))
	{
	    elements = split_line(line,'\t');
	    if (elements.size() != 4)
	    {
		throw std::invalid_argument( "Invalid Map File" );
	    }

	    temp = elements[0];
	    chrom = std::atoi(temp.c_str()) - 1;
	    temp = elements[3];
	    pos = std::atoi(temp.c_str());
	    positions[chrom].push_back(pos);
       }
    }
    else
	throw "Error Reading Map File !";
    return (positions);
}

int range_search(int val, std::vector<int> positions)
{
    /* range search for val in a vector [4,10,15,20]
    for e.g. if val is 6 then will return 0 (b/w 4,10)
    positions is an integer vector in sorted order
    */
    if (val < positions.front() || val > positions.back())
        return -1; // the value was not in any of the ranges    
    else
    {
	auto it = std::lower_bound(positions.begin(),positions.end(),val); 
	return (it - positions.begin() - 1); // get distance from beginning,and return the index in the vector
    }	    
}

void create_ped(const std::string& gvcf_file, const std::string& vcf_file, const:: std::string& map_file, const std::string &out_file, const std::string &sample_name, const std::string &family_name)
{
    /* The output file */
    std::ofstream output;

    /* Read the vcf */
    auto vcf_data = file_reader(vcf_file,"vcf");
    auto vcf_positions = std::get<0>(vcf_data);
    auto vcf_genotypes = std::get<1>(vcf_data);
      
    /* Read the gvcf */
    auto gvcf_data = file_reader(gvcf_file,"gvcf");
    auto gvcf_positions = std::get<0>(gvcf_data);
    auto gvcf_genotypes = std::get<1>(gvcf_data);

    /* Read the map file */    
    std::vector< std::vector<int> >::const_iterator chrom;
    std::vector<int>::const_iterator col;
    auto map_coordinates = return_map_coordinates(map_file);
    int pos;
    std::string geno;
    
    /* Iterate through map coordinates and write in ped format */
    output.open(out_file);    
    int i = 0;
    int wierd = 0;
    if (output)
    {
        output << family_name << " " << sample_name <<  " 0 0 1 0 "; // Add random values for phenotype columns
        for (const auto& chrom : map_coordinates)
	{
	    for (const auto& item : chrom)
	    {
		// Do a range search   
		pos = range_search(item,vcf_positions[i]); // First search the vcf
		if (pos != -1)
		{
		    geno = vcf_genotypes[i][pos];
		}
		else // Try the gvcf
		{
		    pos = range_search(item,gvcf_positions[i]);
		    if (pos != -1)
			geno = gvcf_genotypes[i][pos];
		    else // The variant was not found in the gvcf too !
		    {
			std::cout << i << " " << item << "\n";
			geno = "0 0";
			++wierd;
		    }
		}    
		output << geno << " ";
	    }
	    ++i;
	}
	output << '\n';
   }
    else
	throw "Error Writing to output file !";
    
    output.close();
}

int main(int argc, char* argv[])
{
    if (argc != 7)
    {
      std::cerr << "Usage : " << argv[0] << " <gvcf file> <vcf file> <map file> <output file> <sample_name> <family_name>" << '\n';
    } 
       
    create_ped(argv[1],argv[2],argv[3],argv[4],argv[5],argv[6]);
    return 0;
}

