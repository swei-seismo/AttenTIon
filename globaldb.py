## READ ANTELOPE DATABASE, SET GLOBAL VARIABLES
import os
import sys
sys.path.append( os.environ['ANTELOPE'] + '/data/python' )
from antelope.datascope import *

def dbread(subdb):
    global figdir,dbname,db,dbja,dbass,dbwf,dbresp,dbor,sacdir
    figdir='./specfig_sub'+subdb
    sacdir='/Volumes/igppwei/sacfl/'
##    if subdb=='5':
##        dbname='/P/sopac/mydatabases/lau-alldb/lau-alldb'
##    elif subdb=='2':
##        dbname='/P/caichen/research/Lau/lau_pick/lau-all-GSN-2_done/lau-all-GSN-2'
##    elif subdb=='3' or subdb=='33':
##        dbname='/P/weisq/EQLau/database3/lau_working'
##    elif subdb=='4' or subdb=='44':
##        dbname='/P/weisq/EQLau/database4/lau-all-GSN-4'
##    elif subdb=='1':
##        dbname='/P/sopac/lau1backups/2July13/lau-all-GSN-1'
##    elif subdb=='11':
##        dbname='/P/aadams/mydatabases/intermediate/sub01/lau-all-GSN-1'
##    elif subdb=='22':
##        dbname='/P/aadams/mydatabases/intermediate/sub02/lau-all-GSN-2'
    dbdir='/Volumes/igppwei/LAU2010/picked-subsets/'
    if subdb=='1':
        dbname=dbdir+'lau-all-GSNdb-subset-01'
    elif subdb=='11':
        dbname=dbdir+'lau-all-GSNdb-subset-11'
    elif subdb=='2':
        dbname=dbdir+'lau-all-GSNdb-subset-02'
    elif subdb=='22':
        dbname=dbdir+'lau-all-GSNdb-subset-22'
    elif subdb=='3' or subdb=='33':
        dbname=dbdir+'lau-all-GSNdb-subset-3'
    elif subdb=='4' or subdb=='44':
        dbname=dbdir+'lau-all-GSNdb-subset-4'
    db=dbopen(dbname,'r')
    dbar=db.lookup(table='arrival')
    dbsu=dbar.subset('auth!~/orbassoc_l/')
    dbas=db.lookup(table='assoc')
    dbor=db.lookup(table='origin')
    dbwf=db.lookup(table='wfdisc')
    dbresp=dbwf.lookup(table='sensor')
    dbinst=dbwf.lookup(table='instrument')
    dbsite=dbwf.lookup(table='site')
    dbresp=dbresp.join(dbinst)
    dbresp=dbresp.join(dbsite)
    dbja=dbsu.join(dbas)
    dbass=dbja.join(dbor)
##    if subdb=='22':
##        db2=dbopen(dbdir+'lau-all-GSNdb-subset-02','r')
##        dbwf2=db2.lookup(table='wfdisc')
##        dbresp=dbwf2.lookup(table='sensor')
##        dbinst2=dbwf2.lookup(table='instrument')
##        dbsite2=dbwf2.lookup(table='site')
##        dbresp=dbresp.join(dbinst2)
##        dbresp=dbresp.join(dbsite2)
##        
