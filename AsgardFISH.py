import pandas as pd
import numpy as np

#root file path
root="/Users/ashley/Documents/Research/Post Doctoral Research/AsgardFISHfiles"

#functions

###TestProbe functions###

def tax_filter(df,taxStr,outName,retDf=False,exp=True,rev=False):

    df["taxonomy"]=df.path+df.organismName #create full taxonomy by combining path and organism name

    notTaxStr_indlist=[x for x in df.index if taxStr not in df.taxonomy[x]] #get the numeric indices of all rows where the search string isnt in the taxonomy column
    
    if rev: #if you want to keep everything that is filter term
        TaxStrDf=df.drop(df.index[notTaxStr_indlist])#define new "not taxon" dataframe as all rows from prior searc
    else:   #if you want to keep everything that is not filter term
        TaxStrDf=df.iloc[notTaxStr_indlist] #define new "not taxon" dataframe as all rows from prior search
    
    if exp:
        TaxStrDf.to_csv(root+"/TestProbe/processed/"+outName+".csv",index=False) #export the dataframe to csv file, and set index to false because numeric
    
    if retDf: #if you want the dataframe back in addition to exporting
        return TaxStrDf

def competitor_block(seqStr,pos,base):

    #assumes string structure as "AGGACGGU-=================U===-GGCGCCUCA"
   # print("seqStr",seqStr)
    startSub=seqStr.find("-")+1 #beginning of base position count. cannot use "=" in case base replaces first = from left or right
    endSub=seqStr.rfind("-")-1 #end of base position count
    
    subStr=seqStr[startSub:endSub+1] #define substring as portion only made of bases and "="

    revSubStr=subStr[::-1] #reverse because probe and sequences are reverse compliment style
   # print("pos","reverse substring",pos,revSubStr,"\n")
        
    boolEval=revSubStr[pos-1].upper()==base.upper() #boolean evaluation (true/false). -1 accounts for indexing start at 0
    #I am converting everything to upper case here incase there is a repeat region (indicated by a lowercase/soft masking). Although it signfies a collapsed region, that base is still occurring there.
    
    return boolEval

def comp_diff_binds(probe,compL,probename,exp=True): #draft function to determine inputs for competitor block- requires a probe and a list of competitor probes
    #initialize probename to "" in the case no export is wanted
    
    diffBasesBind=[]#empty list to append base position, base tuples for downstream code processing
    
    if exp:
        basesBindByProbeExp=[] #empty list for exporting the table of competitor and what bases are blocked per competitor
    
    for comp in compL: #for each competitor probe
        if exp:
            diffBasesBindPerProbeExp=[] #list of base/position tuples that will go into the export list by probe
            
        print ("\ncompetitor is", comp)
                
        for count, b in enumerate(probe): #iterate over bases in the probe
            #print("count",count,"probe",b,"comp",comp[count])
        
            if b !=comp[count]: #check same position in competitor probe, if it isn't the same, find what base the competitor is binding to
                #print("\ndoes not match\n")

                if comp[count]=="A":
                    diffBasesBind.append((count+1,"T"))
                    diffBasesBind.append((count+1,"U"))
                    
                    if exp:
                        diffBasesBindPerProbeExp.append((count+1,"T"))
                        diffBasesBindPerProbeExp.append((count+1,"U"))
                    
                elif comp[count]=="T" or comp[count]=="U":
                    diffBasesBind.append((count+1,"A"))
                    
                    if exp:
                        diffBasesBindPerProbeExp.append((count+1,"A"))
                    
                elif comp[count]=="C":
                    diffBasesBind.append((count+1,"G"))
                    
                    if exp:
                        diffBasesBindPerProbeExp.append((count+1,"G"))
                
                elif comp[count]=="G":
                    diffBasesBind.append((count+1,"C"))
                    
                    if exp:
                        diffBasesBindPerProbeExp.append((count+1,"C"))
                    
        if exp:
            str_diffBasesBindPerProbeExp="; ".join(map(str, diffBasesBindPerProbeExp)) #convert the list of tuples to a string to have one entry per probe
            str_diffBasesBindPerProbeExp=str_diffBasesBindPerProbeExp.replace(",",":") #replace the commas with a different delimiters so it doesnt interfere with csv loading
        
            print("blocks for probe",str_diffBasesBindPerProbeExp,"\n")
            basesBindByProbeExp.append((comp,str_diffBasesBindPerProbeExp))
      
    if exp:
        basesBindByProbeExpDict=dict(basesBindByProbeExp) #for exporting. key is probe, values are list of base/position tuples
        print("collated for export",basesBindByProbeExpDict,"\n")
    
        basesBindByProbeExpDf=pd.DataFrame(basesBindByProbeExpDict.items(),columns=["Competitor","position: base blocked"])
        basesBindByProbeExpDf.set_index("Competitor",inplace=True)
    
        print("df",basesBindByProbeExpDf,"\n")

        basesBindByProbeExpDf.to_csv(root+"/TestProbe/processed/probecompetitors_baseblocks_"+probename+".csv",index=True)
    
    return diffBasesBind #returns a list of position, bound-base tuples. For downstream code processing.

def apply_comptupes_to_hl(hl,compTupLst,taxStr,probename,exp=True, rev=False,hits=False,ret=True): #enter dataframe of hitlist/sequences and the list of tuples with base position/what base competitor binds to

    if rev:
        revStr="BoundBy"
    else:
        revStr="After"
        
    if hits:
        hitsStr=""
    else:
        hitsStr="not"

    for count,t in enumerate(compTupLst): #each iteration through the position/bound base tuple list adds a boolean column (is the binding base blocked or not)
        print(t,count)
        hl["compBlock"+str(count)]=hl.sequenceMatch.apply(lambda x: competitor_block(x,t[0],t[1]))

    allBoolCols=[i for i in hl.columns if "compBlock" in i] #locate all boolean analysis columns
    hl["allBools"]=np.sum(hl[allBoolCols],axis=1)#sum all boolean columns. 1= True+False or True+True. 0=False+False. We want to drop everything but where it is 0
    
    print("boolean sums\n",hl["allBools"])
    
    if rev:
        hlAfterComp=hl.loc[hl.allBools!=0] #filter dataframe accordingly if you want to export what sequences the competitor does bind to
    else:
        hlAfterComp=hl.loc[hl.allBools==0] #filter dataframe accordingly
       
    print(hlAfterComp)
    dropCols=[i for i in hlAfterComp.columns if "compBlock" in i or i=="allBools"] #Drop all boolean columns before exporting

    hlAfterComp.drop(dropCols, axis=1,inplace=True)
    hl.drop(dropCols, axis=1,inplace=True)
    
    if exp:
        hlAfterComp.to_csv(root+"/TestProbe/processed/"+hitsStr+taxStr+"_"+probename+"_1mm_"+revStr+"Com.csv",index=False)
    
    if ret:
        return hlAfterComp

def unique_accesions_TestProbe(testProbeDf,outName="",extPath="",exp=True,retList=False): #exports and/or returns a list of unique accessions from a TestProbe dataframe

    accessions=testProbeDf.primaryAccession.tolist()
    uniqueAccessions=np.unique(accessions).tolist()
            
    if exp:
        writeListToFile(uniqueAccessions,outName,extPath)

    if retList:
        return uniqueAccessions


###FASTA file functions###

def read_fasta_to_list(fileName): #reads in fasta file as one large string, converts to entry list

     fastaFile=open(fileName) #open file
     fastaAsString=fastaFile.read() #read in text as one large string
     fastaFile.close()
     
     fastaAsList=fastaAsString.split(">") #split string into a list using fasta header as delimiter
     del fastaAsList[0] #delete first item which will just be a space
     
     return fastaAsList
     

def unique_accesions_fasta(fastaList,outName="",exp=True,retList=False): #creates a list of unique accessions from fasta file. Can return file as a list and/or export list to text file

    allAccessions=[]
    
    for item in fastaList: #in each fasta entry
        entryAsList=item.split(" ") #split the entry into a list using space as delimiter
        accession=entryAsList[0] #text string before first space is the accession
        print(accession)
        allAccessions.append(accession) #append accession to list
        
    uniqueAccessions=np.unique(allAccessions).tolist() #cut down to unique accessions only
        
    if exp:
        writeListToFile(uniqueAccessions,outName)


    if retList:
        return uniqueAccessions

def writeListToFile(outList,outName,extPath="/SILVA138_refNR/out/",nl="\n",fileExt="txt"): #write list to text file. Leave new line character as an option so if exporting list entries with \n already in, as with fasta, can bypass

    outFile = open(root+extPath+outName+"."+fileExt, "w") #create outfile
    
    for element in outList:
        outFile.write(element + nl)
        
    outFile.close()

def read_linedelimited_text_to_list(filePath):
    delimitedFile = open(filePath+".txt", "r")
    delimitedFileAsList = delimitedFile.readlines()
    
    delimitedFileAsList=[x.replace("\n","") for x in delimitedFileAsList]
    
    return delimitedFileAsList
    
def filter_fasta_by_accessions(accessionL,fastaAsList,outname,extPath="/fasta_targetAccessions/"): # given a list of accessions and a fasta file, export a fasta containing only the accessions in the input list

    matchedEntries=[]
    checkL=[]

    for a in accessionL: #per accession
        #print("accession",a)
        for entry in fastaAsList: #iterate through fasta file entries
            list_entry=entry.split(" ") #split an entry into a list using space as delimiter
            #print("accession as list",list_entry)
            longAccession=list_entry[0] #define the full length accession number as the first item before a space
            #print("long accession number",longAccession)
            list_longAccession=longAccession.split(".") #split the full length accession into components using "." as the delimiter
            #print("long accession number as list",list_longAccession)
            primaryAcc=list_longAccession[0] #primary accession is characters before the first period. This is required because test probe accessions are primary accession format.
            #print("primary accession",primaryAcc)
        
            if a==primaryAcc: #if accessions match
                #print("accessions match",a,primaryAcc)
                matchedEntries.append(">"+entry)
                checkL.append(primaryAcc)
    
            #print("\n\n")

        #print("\n\n\n")

    writeListToFile(matchedEntries,outname,extPath+"fastas/",nl="",fileExt="fasta") #output fasta file
    writeListToFile(checkL,outname+"_checkListAccessions",extPath+"checkfiles/") #output the list of matched accessions as a test file to check work
    
    
###Runners###
def run_fasta_target_analysis(taxStr,fastaFileName,retAcc=True,retFastaAsList=False): #required to run before testProbe for coverage calculations
    fastaAsList=read_fasta_to_list(fastaFileName) #read in a fasta file and turn it into a list of accessions
    uniqueAccessions=unique_accesions_fasta(fastaAsList,outName=taxStr+"UniqueAccessions",exp=True,retList=True) #export and return a list of unique accessions from the fasta file
    
    if retAcc and not retFastaAsList:
        return uniqueAccessions
    
    elif retAcc and retFastaAsList: #retFastaAsList also returns the silva fasta file as a list of entries
        return uniqueAccessions,fastaAsList
        
    elif not retAcc and retFastaAsList:
        return fastaAsList
    
def run_full_testprobe_analysis_phase1(hitdf,taxStr,probename,probe,compL=[],retHLs=True,comps=True,retTargetAcc=False):
    #isolate non-targets
    notTax=tax_filter(hitdf,taxStr,"not"+taxStr+"_"+probename+"_1mm",retDf=True) #example , not targets filter
    Tax=tax_filter(hitdf,taxStr,taxStr+"_"+probename+"_1mm",rev=True,retDf=True) #example , targets filter
    #print("targets hitlist",Tax,"\n")

    if comps:
        #competitor analysis
        #loop through using a list of position, base tuples...position indicates position from the left of the probe sequence and base indicates the base that is bound by the probe
        compBindTupList=comp_diff_binds(probe,compL,probename)
        print("compBindTupList",compBindTupList,"\n")

        #iterate through position/bases bound and a boolean column (is the binding base blocked or not), keep whatever isn't blocked from target and non-target dataframes
        notTaxAfterComp=apply_comptupes_to_hl(notTax,compBindTupList,taxStr,probename) #non-targets
        TaxAfterComp=apply_comptupes_to_hl(Tax,compBindTupList,taxStr,probename,hits=True) #Targets
        
        #export unique list of target accessions that bind to the probe after competitors to pull from SILVA fasta files for alignments
        unique_accesions_TestProbe(TaxAfterComp,outName="accessions"+"_"+taxStr+"_"+probename+"_1mm_AfterCom",extPath="/TestProbe/processed/",exp=True,retList=False)
    
        #export target and non-target sequences that the competitor(s) bind(s) to
        apply_comptupes_to_hl(notTax,compBindTupList,taxStr,probename,rev=True,ret=False)
        apply_comptupes_to_hl(Tax,compBindTupList,taxStr,probename,rev=True,hits=True,ret=False)
    
        if retHLs and not retTargetAcc:
            return notTax,Tax,notTaxAfterComp,TaxAfterComp
            
        elif retHLs and retTargetAcc:
            return unique_accesions_TestProbe,notTax,Tax,notTaxAfterComp,TaxAfterComp
        
    else:
        #export unique list of target accessions that bind to the probe to pull from SILVA fasta files for alignments
        unique_accesions_TestProbe(Tax,outName="accessions"+"_"+taxStr+"_"+probename+"_1mm",extPath="/TestProbe/processed/",exp=True,retList=False)
        
        if retHLs and not retTargetAcc:
            return notTax,Tax
            
        elif retHLs and retTargetAcc:
            return unique_accesions_TestProbe, notTax,Tax
    
def calculate_coverage_specificity(targetsPreComp,nonTargetsPreComp,uniqueAccessions,taxStr,probename,targetsPostComp="",nonTargetsPostComp="",exp=False,ret=True):

    preCompCoverage=len(targetsPreComp)/len(uniqueAccessions)*100 #coverage is number of targets hit by probe relative to total unique accessions in database
    specificityPreComp=len(targetsPreComp)/(len(nonTargetsPreComp)+len(targetsPreComp))*100 #specificity is number of targets hit by probe relative to total targets+non targets
    
    if isinstance(targetsPostComp,pd.DataFrame): #if there are competitors, targetsPostComp will be a dataframe rather than being initialized as an empty string
        postCompCoverage=len(targetsPostComp)/len(uniqueAccessions)*100
        specificityPostComp=len(targetsPostComp)/(len(nonTargetsPostComp)+len(targetsPostComp))*100
    else:
        postCompCoverage=np.nan
        specificityPostComp=np.nan
        
    specCovDf=pd.DataFrame([(probename,taxStr,preCompCoverage, postCompCoverage, specificityPreComp, specificityPostComp)],columns=["probe","target","pre-competitor coverage","post-competitor coverage","pre-competitor specificity","post-competitor specificity"])
    
    if exp:
        specCovDf.to_csv(root+"/specificity_coverage/"+taxStr+"_"+probename+"speccov_1mm.csv",index=False)
    if ret:
        return specCovDf

def create_subfasta_from_proc_target_accessions(taxStr,probename,targetFastaAsList,compStr="_AfterCom"): #re-opens text file of target accessions bound by probe (after competitors if applicable) and generates fasta file of those accessions

    taxPostCompAcc=read_linedelimited_text_to_list(root+"/TestProbe/processed/accessions_"+taxStr+"_"+probename+"_1mm"+compStr)
    filter_fasta_by_accessions(taxPostCompAcc,targetFastaAsList,taxStr+"_"+probename+"_1mm"+compStr+"Fasta",extPath="/fasta_targetAccessions/")
    
###body###

##data set one- all Loki

LokiSILVA138Fasta=root+"/SILVA138_refNR/searchresults/arb-silva.de_2021-10-27_id1073105_withoutgaps_Lokiarchaeia.fasta"
Loki1378_1mm=pd.read_csv(root+"/TestProbe/searchresults/arb-silva.de_testprobe_hitlist_1048769_LOK1378_1mm.csv",delimiter=";")
Loki1378Pr="GTGTGCAAGGAGCAGAGACG"
c16Loki1378="GTGTGCAAGGAGCAGGGACG" #competitor
c20Loki1378="GTGTGCAAGGAGCAGGGACA" #competitor
c2Loki1378="GAGTGCAAGGAGCAGGGACG" #competitor

lokiarchaeiaUniqueAcc,lokiarchaeiaFastaAsList=run_fasta_target_analysis("Lokiarchaeia",LokiSILVA138Fasta,retAcc=True,retFastaAsList=True)
notLokiarchaeiaPreComp,lokiarchaeaPreComp,notLokiarchaeiaPostComp,lokiarchaeiaPostComp=run_full_testprobe_analysis_phase1(Loki1378_1mm,"Lokiarchaeia","Loki1378",Loki1378Pr,[c16Loki1378,c20Loki1378,c2Loki1378])
calculate_coverage_specificity(lokiarchaeaPreComp,notLokiarchaeiaPreComp,lokiarchaeiaUniqueAcc,"Lokiarchaeia","Loki1378",lokiarchaeiaPostComp,notLokiarchaeiaPostComp,exp=True,ret=False)
create_subfasta_from_proc_target_accessions("Lokiarchaeia","Loki1378",lokiarchaeiaFastaAsList)

##data set two- Loki clade A
##
Loki1183_1mm=pd.read_csv(root+"/TestProbe/searchresults/arb-silva.de_testprobe_hitlist_1048765_LOK1183_1mm.csv",delimiter=";")
Loki1183Pr="CTGACCTGCCTTTGCCCGCT"
Loki1183Cp="CTGACCTGCCGTTGCCCGCT" #competitor

notLokiarchaeiaAPreComp,lokiarchaeaAPreComp,notLokiarchaeiaAPostComp,lokiarchaeiaAPostComp=run_full_testprobe_analysis_phase1(Loki1183_1mm,"Lokiarchaeia","Loki1183",Loki1183Pr,[Loki1183Cp])
calculate_coverage_specificity(lokiarchaeaAPreComp,notLokiarchaeiaAPreComp,lokiarchaeiaUniqueAcc,"Lokiarchaeia_A","Loki1183",lokiarchaeiaAPostComp,notLokiarchaeiaAPostComp,exp=True,ret=False)
create_subfasta_from_proc_target_accessions("Lokiarchaeia","Loki1183",lokiarchaeiaFastaAsList)

#data set three- Heimdallarchaeia
HeimdalSILVA138Fasta=root+"/SILVA138_refNR/searchresults/arb-silva.de_2021-12-07_id1092153_withoutgaps_Heimdallarchaeia.fasta"
Heim329_1mm=pd.read_csv(root+"/TestProbe/searchresults/arb-silva.de_testprobe_hitlist_1048771_ HEIM329_1mm.csv",delimiter=";")
Heim329Pr="GCACTCGCAGAGCTGGTTTTACC"

heimdallarchaeiaUniqueAcc,heimdallarchaeiaFastaAsList=run_fasta_target_analysis("Heimdallarchaeia",HeimdalSILVA138Fasta,retAcc=True,retFastaAsList=True)
notHeimdallarchaeiaPreComp,HeimdallarchaeiaPreComp=run_full_testprobe_analysis_phase1(Heim329_1mm,"Heimdallarchaeia","Heim329",Heim329Pr,comps=False)
calculate_coverage_specificity(HeimdallarchaeiaPreComp,notHeimdallarchaeiaPreComp,heimdallarchaeiaUniqueAcc,"Heimdallarchaeia","Heim329",exp=True,ret=False)
create_subfasta_from_proc_target_accessions("Heimdallarchaeia","Heim329",heimdallarchaeiaFastaAsList,compStr="")
####

