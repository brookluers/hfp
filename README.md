--- basic.go ---
fits the proportional hazards models



--- data.go ---
reads hfdat.gob.gz and creates binary column files

Elix_0, Elix_1... are 0/1 indicator variables for each of the Elixhauser categories

TG_01, ... drug theraputic group values

time, DOB, gender, etc...

--- hfdat.go --- 
Read from Marketscan, convert go Gob

read the A, O, S, I, F, D files from MarketScan, segmenting by Enrolid (rows sorted by Enrolid then date)
     into dstreams

concurrently process each bucket of data

dobucket(k int) processes the kth bucket
	   joins the dstreams returned by setupBucket
	   
	   loops over the subjects in each bucket
	   
	   checking for eigibility, compiling drug/procedures, etc
	   checking if heart failure present (in obs. period or in prediction period)
	   
	   keep all cases (heart failure present in prediction period)
	   keeping random 10% sample of controls
	   	   
	   stores each retained subject in a utils.Drec struct
	   passes the Drec into the rslt channel
	   
setupBucket(k int) returns an array of dstreams, one dstream for each A, I, O...


harvest()
	ranges through the rslt channel
	writes all the records to hfdat.gob
	


--- reduce.go ---
extracts 20 factors from the procedure codes using SVD

--- kshedden/gocols/config repository ---
parse configuration of column-stored compressed data 

e.g. for a given bucket number, return the path to the data in that bucket
     return the codes for a factor variable (map from strings to integers)
     read a json of data types for each column



--- utils/defs.go --- 
Drec struc represents a single person
     indicator of heart failure
     sex
     elixhauser indicators
     drug group indicators
     procedure codes
