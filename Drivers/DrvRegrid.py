#!/usr/bin/env python
# Import packages 
import sys
import argparse as arg

import time
import numpy as np

# import modules in other directories
sys.path.append('../Regridder/')
import GenRegrid as GnR
import Initialize as Init
import WriteDST as Wrt
import ReadInSrc as Rd 

import importlib
importlib.reload( Rd )

#Rdry = Con.Rdry() # 



def main(year, month, day, hour, Dst, DstVgrid, Src, IC_for_pg, RegridMethod = 'CONSERVE'):
    import calendar
    
    tic_total = time.perf_counter()
    days_in_month = calendar.monthrange(year,month)[1]

    print( f"About to process {year:n}-{month:n}-{day:n}")

    #RegridMethod = 'CONSERVE'
    lnPS=False
    if(lnPS==True):
        ver='lnPS'
    else:
        ver=''

    # Override this when it can't be true
    if (Dst not in ('ne480np4','ne240np4','ne120np4','ne30np4')):
        IC_for_pg = False
        print(f'  Setting IC_for_pg={IC_for_pg} because cannot be otherwise for {Dst} ' )
    else: 
        print(f'  {IC_for_pg}: This is making an IC file for {Dst} ' )


    ret1 = Init.prep(Dst=Dst, DstVgrid=DstVgrid ,Src=Src, WOsrf=True, RegridMethod=RegridMethod , IC_for_pg=IC_for_pg )
    sys.stdout.flush()
    if (day==99):
        for iday in np.arange( days_in_month):
            ret2 = Rd.get_Src( year=year ,month=month ,day=iday+1 , hour0=99 )
            sys.stdout.flush()
            ret3 = GnR.xRegrid(HorzInterpLnPs=lnPS )
            sys.stdout.flush()
            ret4 = Wrt.write_netcdf(version=ver) #+'Test01')

    else:
        ret2 = Rd.get_Src( year=year ,month=month ,day=day , hour0=hour )
        sys.stdout.flush()
        ret3 = GnR.xRegrid(HorzInterpLnPs=lnPS )
        sys.stdout.flush()
        ret4 = Wrt.write_netcdf(version=ver) #+'Test01')
        
    code = 1
    toc_total = time.perf_counter()

    pTime = f"Total processing time was  {toc_total - tic_total:0.4f} seconds"
    print(pTime)

if __name__ == "__main__":
    # argument: indir -> get all nc files in this directory
    # argument: map -> the offlinemap file already prepared
    # argument: outdir -> directory where remapped files should go
    # my_parser = arg.ArgumentParser()
    # my_parser.add_argument("--month", type=int)
    # my_parser.add_argument("--year", type=int)
    # args = my_parser.parse_args()
    
    my_parser = arg.ArgumentParser()
    my_parser.add_argument("--month",    type=int, default=1)
    my_parser.add_argument("--year",     type=int, default=2010)
    my_parser.add_argument("--day",      type=int, default=1)
    my_parser.add_argument("--hour",     type=int, default=99)
    my_parser.add_argument("--Dst",      type=str, default="ne30pg3")
    my_parser.add_argument("--Src",      type=str, default="ERA5")
    my_parser.add_argument("--DstVgrid", type=str, default="L58")
    args = my_parser.parse_args()
    main(args.year, args.month, args.day, args.hour, args.Dst, args.DstVgrid, args.Src )
