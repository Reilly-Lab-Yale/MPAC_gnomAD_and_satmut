#!/usr/bin/env python
# coding: utf-8

# This notebook annotates all gnomad variants with their corresponding malinouis predictions, the phyloP scores associated with their genomic locations, and the genomic regions (enhancer or not) they fall within. 
# 
# (Execute notebook after crunching wig into csv, as per other file)
# 
# ## Import relevant libraries

# In[1]:


from pyspark import SparkConf, SparkContext
from pyspark.sql import SparkSession
from pyspark.sql.types import StructType, StructField, StringType, IntegerType, FloatType
import pyspark.sql.functions as F
import os


# ## Create spark session

# In[2]:


if 'spark' in locals() and spark!=None:
    spark.stop()

    #are we running the actual script, or just testing?
for_real=True

spark=None

spark = SparkSession.builder.getOrCreate()


# ## select chromosome

# In[18]:


chromosome="None"

if "which_chr" in os.environ:
    chromosome=os.environ['which_chr']
    print("only crunching chromosome "+chromosome)
else:
    print("error : did not find which chromosome we are supposed to crunch")
    exit(-1)


# ## load in data

# In[19]:


#define the phylop tsv schema
phylop_schema = StructType([
    StructField("CHROM", StringType(), True),
    StructField("POS", IntegerType(), True),
    StructField("P_ANNO", FloatType(), True),
])

#read in the phylop tsv
phylop_anno = spark.read \
    .option("comment", "#") \
    .option("delimiter", "\t") \
    .schema(phylop_schema) \
    .csv("/home/mcn26/varef/scripts/noon_data/1.0.format_zoonomia_phylop/out_processed.tsv", header=False)


# In[21]:


vcf_schema = StructType([
    StructField("CHROM", StringType(), True),
    StructField("POS", IntegerType(), True),
    StructField("ID", StringType(), True),
    StructField("REF", StringType(), True),
    StructField("ALT", StringType(), True),
    StructField("QUAL", StringType(), True),
    StructField("FILTER", StringType(), True),
    StructField("INFO", StringType(), True),

])

vcf = spark.read \
    .option("comment", "#") \
    .option("delimiter", "\t") \
    .schema(vcf_schema) \
    .csv(f"/home/mcn26/varef/scripts/noon_data/0.merge/combined.{chromosome}.vcf.gz", header=False)


# In[23]:


#load the encode enhancer dataset

#define genomic region annotation (bed file) schema
bed_schema = StructType([
    StructField("CHROM", StringType(), True),
    StructField("START", IntegerType(), True),
    StructField("STOP", StringType(), True),
    StructField("misc_1", StringType(), True),
    StructField("misc_2", StringType(), True),
    StructField("REGION_TYPE", StringType(), True),
])


all_CREs=spark.read \
    .schema(bed_schema) \
    .option("delimiter", "\t") \
    .csv("/home/mcn26/varef/data/ENCODE/SCREEN_v4_cCREs_agnostic/GRCh38-cCREs.V4.bed.gz")

all_CREs = all_CREs.drop("misc_1", "misc_2")


# ## Split relevant information from the VCF info field
# ...into its own columns.

# In[4]:


####The `INFO` field contains a lot of useful information, but it is all smashed together into a string. 
#Let's extract information from that string. 

keys_to_extract = [#NONE CAN BE SUBSTRINGS OF THE OTHERS
    "K562__ref", "HepG2__ref", "SKNSH__ref", "K562__alt", "HepG2__alt", "SKNSH__alt",
    "K562__skew", "HepG2__skew", "SKNSH__skew", "AC", "AN", "AF", "cadd_phred",# "P_ANNO" already in its own column
]

# Apply the regexp_extract function to the DataFrame to create new columns for each key.
# The expression '([^;]*)' captures any sequence of characters that are not a semicolon,
# which is assumed to be the delimiter for the key-value pairs in the 'INFO' column.

for key in keys_to_extract:

    #df = df.withColumn(key, regexp_extract(col("INFO"), "{}=([^;]+);?".format(key), 1))
    #when we find something put it, whne we don't put None
    vcf = vcf.withColumn(key, 
                       F.when(
                           F.regexp_extract(F.col("INFO"), "{}=([^;]+);?".format(key), 1) != "",
                           F.regexp_extract(F.col("INFO"), "{}=([^;]+);?".format(key), 1)).otherwise(None))


# ## Filter out the gorp
# 
# - Filter out the gorp
#     - Remove low-quality variants that don't pass GNOMAD's own filters.
#     - Remove low-quality variants not queried in a large number of individuals
#     - Remove variants with a MAF of 0 (they don't really "vary") in the population if they dont exist.
#     - Remove those variants that don't have all of the relevant metrics.
#     - Remove non-SNP variants

# In[5]:


vcf = vcf.filter(
    #make sure we have the necessary population stats
    (F.col("AF").isNotNull()) & 
    (F.col("AC").isNotNull()) & 
    (F.col("AN").isNotNull()) & 
    
    #do NOT make sure we a cadd_phred score : we don't care
    #(F.col("cadd_phred").isNotNull()) & 

    #check variant has been queried in a reasonably large number of people
    #same as gnomad's own warning threshold
    #see readme
    (F.col("AN").cast("int") >= 76156) &
    
    #gnomad filters passed. See original gnomad vcf header for spec.
    (F.col("FILTER") == "PASS") 
    
)


# ## Add the genomic region annotations

# In[15]:


def add_genomic_annotation(loci,regions,name):
    

    # as usual, we have to worry about coordinate systems
    # VCFs are 1-based
    # BEDs are 0-based
    # Chr1        T   A   C   G   T
    #           | | | | | | | | | |
    # 1 based   | 1 | 2 | 3 | 4 | 5
    # 0 based   0   1   2   3   4

    #subset the provided regions to only those of the region of interest
    regions_subset = regions.filter( (F.col("REGION_TYPE") == name ) )
    
    #joins 
    result = loci.join(
        regions_subset,
        (loci.CHROM == regions_subset.CHROM) & 
        (loci.POS > regions_subset.START) & 
        (loci.POS <= regions_subset.STOP),
        how='left'
    )
    
    # note the potential for a bug here. Left join can perform one-to-many
    # and will create multiple entries. That is, if a variant falls inside
    # multiple overlapping regions of the same region type, the variant will be duplicated. 
    # However, as I've checked in helper/check_overlap/, this should not be a problem
    # for this particular dataset
    
    # Add a boolean column 'is_in_region'
    result = result.withColumn("is_in_"+name, F.col("start").isNotNull())

    # Select only columns from df_loci and the new boolean column
    final_result = result.select(loci["*"], "is_in_"+name)
    
    return final_result


# In[7]:


#get all the unique values we might want to add
unique_values=all_CREs.select("REGION_TYPE").distinct()

unique_list = [row['REGION_TYPE'] for row in unique_values.collect()]


# In[18]:


#make a boolean column for each genomic region annotation
genomic_regions_added=vcf
for region in unique_list:
    genomic_regions_added=add_genomic_annotation(loci=genomic_regions_added,regions=all_CREs,name=region)


# ## Add the phyloP annotations

# In[31]:


phyloP_annotated = genomic_regions_added.join(phylop_anno, on=["CHROM", "POS"], how="inner")
#inner join means we'll drop all those variants without scores...


# Compute mean malinois reference activity and mean malinois skew

# In[32]:


#reference activity

phyloP_annotated = phyloP_annotated.withColumn("K562__ref", F.col("K562__ref").cast("float"))
phyloP_annotated = phyloP_annotated.withColumn("HepG2__ref", F.col("HepG2__ref").cast("float"))
phyloP_annotated = phyloP_annotated.withColumn("SKNSH__ref", F.col("SKNSH__ref").cast("float"))

phyloP_annotated=phyloP_annotated.withColumn("mean_ref", (F.col("K562__ref") + F.col("HepG2__ref") + F.col("SKNSH__ref")) / 3)

#skew
phyloP_annotated = phyloP_annotated.withColumn("K562__skew", F.col("K562__skew").cast("float"))
phyloP_annotated = phyloP_annotated.withColumn("HepG2__skew", F.col("HepG2__skew").cast("float"))
phyloP_annotated = phyloP_annotated.withColumn("SKNSH__skew", F.col("SKNSH__skew").cast("float"))

phyloP_annotated=phyloP_annotated.withColumn("mean_skew", (F.col("K562__skew") + F.col("HepG2__skew") + F.col("SKNSH__skew")) / 3)


# ## Add rarity codes

# In[33]:


# First, convert 'AC' and 'AF' to integer and float types respectively using withColumn
df = phyloP_annotated.withColumn("AC", F.col("AC").cast("int"))
df = df.withColumn("AF", F.col("AF").cast("float"))

# Calculate MAF using the Spark SQL functions
df = df.withColumn("MAF", F.least(F.col("AF"), 1 - F.col("AF")))


df = df.withColumn("category", 
                   F.when((F.col("MAF") == 0.0) | (F.col("AC") == 0), "MAF_OR_AC_IS_ZERO")
                    .when(F.col("AC") == 1, "SINGLETON")
                    .when(F.col("MAF") < 0.01 / 100, "ULTRARARE")
                    .when(F.col("MAF") < 0.1 / 100, "RARE")
                    .when(F.col("MAF") < 1 / 100, "LOW_FREQ")
                    .otherwise("COMMON"))


# In[34]:


df.write.option("codec", "org.apache.hadoop.io.compress.GzipCodec") \
    .option("delimiter", "\t") \
    .csv(f"/home/mcn26/varef/scripts/noon_data/2.0.annotate/annotated_output_{chromosome}.csv.gz", header=True, mode="overwrite")

