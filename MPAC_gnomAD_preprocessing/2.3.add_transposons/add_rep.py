#!/usr/bin/env python
# coding: utf-8

# In[2]:


from pyspark.sql import SparkSession
from pyspark.sql import functions as F
import pyspark.sql.types as T
import os


# In[3]:


spark = SparkSession.builder.appName("add_TE").getOrCreate()


# In[4]:


chromosome="NONE"
chromosome="chr22"

if "which_chr" in os.environ:
    chromosome = os.environ['which_chr']

if chromosome=="NONE":
    print("error : did not find which chromosome we are supposed to crunch!")
    exit(-1)
else:
    print("only crunching chromosome "+chromosome)


# In[5]:


variant_path=f"/home/mcn26/varef/scripts/noon_data/indel/2.0.annotate/annotated_output_{chromosome}.csv.gz/*.csv.gz"

variants=spark.read.option("delimiter","\t") \
    .csv(variant_path, header=True, inferSchema=True)


# In[6]:


transposon_path="/home/mcn26/varef/data/repmask/hg38.fa.out.fixed.gz"


# This file is a modified version of the output of repeatmasker. It's a very strange format, which I have munged slightly to be more tractable. Take a look at `helper/prechew_transposons` for more information.

# In[7]:


transposon = spark.read \
    .option("sep", "\t") \
    .option("header", "true") \
    .csv(transposon_path)


# In[8]:


transposon.take(3)


# In[24]:


# not terrifically efficient, should probably be pre-processes for
# faster speeds, but this is a quick step, so it shouldn't matter terribly.
transposon_chr=transposon.filter(F.col("query_sequence") == chromosome)


# In[26]:


columns_to_keep = ["query_sequence","pos_in_query_begin","pos_in_query_end"]
transposon_chr = transposon_chr.select(columns_to_keep)


# In[38]:


# the variants is 1-based, as it is from a VCF
# The transposon table is 1-based, 
# as we discovered in /helper/prechew_transposons/1.0.determine_base
# this makes checking ranges quite easy !

#Chr1        T   A   C   G   T
#          | | | | | | | | | |
#1 based   | 1 | 2 | 3 | 4 | 5
#0 based   0   1   2   3   4

#broadcast transposon for performance
transposon_broadcast = F.broadcast(transposon_chr)

#left-join
#note that this will duplicate variants if they fall into multiple transposons
joined=variants.join(transposon_broadcast,
            (variants.POS > transposon_chr.pos_in_query_begin) & 
            (variants.POS <= transposon_chr.pos_in_query_end),
            how='left'
    )

#add a column that records whether each variant is in a transposon or not
result = joined.withColumn("in_rep", 
             F.when(F.col("pos_in_query_begin").isNotNull(), F.lit(True)).otherwise(F.lit(False)))

#drop rep_element columns we don't care about anymore
#important for dedup
result=result.drop(*transposon_broadcast.columns)

# deduplicate variants by grouping by all columns that aren't the newly added in_rep, and taking the max
# max wil thus indicate "in any rep element"
group_columns = [col for col in result.columns if col != "in_rep"]
result = result.groupBy(*group_columns).agg(F.max("in_rep").alias("in_rep"))



# In[39]:


result.write.option("codec", "org.apache.hadoop.io.compress.GzipCodec") \
    .option("delimiter", "\t") \
    .csv(f"/home/mcn26/varef/scripts/noon_data/indel/2.3.add_transposons/{chromosome}.csv.gz", header=True, mode="overwrite")

