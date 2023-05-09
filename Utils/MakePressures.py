import numpy as np
import MyConstants as Con
# Physical Constants
Rgas = Con.Rdry() # 287.0 # J K-1 kg-1
grav = Con.grav() # 9.8

def Pressure ( am, bm, ai, bi, ps , p_00=100_000., Gridkey='tzc' ):
    # Should ASSERT that am ,bm, ai, and bi are 1D

    if ( Gridkey == 'zc' ):
        ncol = np.shape( ps )
        nz = np.shape( am )[0]
        delp = np.zeros( (nz,ncol ) )
        pmid = np.zeros( (nz,ncol ) )
        pint = np.zeros( (nz+1,ncol ) )
        for L in np.arange( nz ):
            pmid[L,:] = am[L]*p_00 + bm[L]*ps[:]
        for L in np.arange( nz+1 ):
            pint[L,:] = ai[L]*p_00 + bi[L]*ps[:]
        for L in np.arange( nz ):
            delp[L,:] = pint[L+1,:]-pint[L,:]
 
    if ( Gridkey == 'tzc' ):
        nt, ncol = np.shape( ps )
        nz = np.shape( am )[0]
        delp = np.zeros( (nt,nz,ncol ) )
        pmid = np.zeros( (nt,nz,ncol ) )
        pint = np.zeros( (nt,nz+1,ncol ) )
        for i in np.arange( nt ):
            for L in np.arange( nz ):
                pmid[i,L,:] = am[L]*p_00 + bm[L]*ps[i,:]
            for L in np.arange( nz+1 ):
                pint[i,L,:] = ai[L]*p_00 + bi[L]*ps[i,:]
            for L in np.arange( nz ):
                delp[i,L,:] = pint[i,L+1,:]-pint[i,L,:]
                        
    if ( Gridkey == 'zyx' ):
        ny, nx = np.shape( ps )
        nz = np.shape( am )[0]
        delp = np.zeros( (nz,ny,nx ) )
        pmid = np.zeros( (nz,ny,nx ) )
        pint = np.zeros( (nz+1,ny,nx ) )
        for L in np.arange( nz ):
            pmid[L,:,:] = am[L]*p_00 + bm[L]*ps[:,:]
        for L in np.arange( nz+1 ):
            pint[L,:,:] = ai[L]*p_00 + bi[L]*ps[:,:]
        for L in np.arange( nz ):
            delp[L,:,:] = pint[L+1,:,:]-pint[L,:,:]
                        
    if ( Gridkey == 'tzyx' ):
        nt, ny, nx = np.shape( ps )
        nz = np.shape( am )[0]
        delp = np.zeros( (nt,nz,ny,nx ) )
        pmid = np.zeros( (nt,nz,ny,nx ) )
        pint = np.zeros( (nt,nz+1,ny,nx ) )
        for i in np.arange( nt ):
            for L in np.arange( nz ):
                pmid[i,L,:,:] = am[L]*p_00 + bm[L]*ps[i,:,:]
            for L in np.arange( nz+1 ):
                pint[i,L,:,:] = ai[L]*p_00 + bi[L]*ps[i,:,:]
            for L in np.arange( nz ):
                delp[i,L,:,:] = pint[i,L+1,:,:]-pint[i,L,:,:]
                        
 
 
    return pmid,pint,delp



def TandP150 ( te, pmid, delp, search_up_L=10, Gridkey='tzc', findHt=150. ):
    # Function finds vertical index of first mid level abov 150m (or findHt)

    # delp is +ve
    # search_up_L was introduced for 'efficiency' - as usual efficiency is a 
    # misguided goal, but here it is ... really complicates code.
    """
    Rough guide to indexing gymnastics necessitated by pompous effort at 'efficiency':
    The argument search_up_L=NS limits search to bottom NS levels of model. In these 
    levels the height above ground level 'zaglm' is calculated.  Then k150x is assigned
    a value of 1 where zaglm<150. This is then summed along Z-axis, which gives an inverted
    index for levels < 150, i.e., the index of the lowest zaglm>150 is:
               L150x = NS - SUM(k150x) - 1
    Then this index in the limited search range has to be converted to an index in the
    full model vertcal grid:
               L150 = nz - 1 + L150x - (NS -1)
    """


    if ( Gridkey == 'zc' ):
        nz, ncol = np.shape( pmid )
        # Set up necessary arrays
        zagl  = np.zeros( (search_up_L+1, ncol) )
        dzagl = np.zeros( (search_up_L, ncol) )
        zaglm = np.zeros( (search_up_L, ncol) )
        k150x = np.zeros( (search_up_L, ncol) , dtype=np.int )
        te150 = np.zeros( (ncol) )
        pm150 = np.zeros( (ncol) )
        for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
            L = nz - 1 + (Lx - (search_up_L-1) )
            dzagl[Lx,:] = ( Rgas *te[L,:]/(grav*pmid[L,:]) ) * delp[L,:]
        for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
            zagl[Lx,:] = zagl[Lx+1,:] + dzagl[Lx,:] 
        for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
            zaglm[Lx,:] = 0.5*(zagl[Lx+1,:] + zagl[Lx,:] )

        # Insane indexing gymnastics for 'efficiency'.
        # First line is straightforward
        k150x = np.where( zaglm < 150. , 1, 0 )
        k150  = search_up_L - np.sum( k150x , axis=0 ) -1
        L150  = nz - 1 + (k150 - (search_up_L-1) )

        # Enforce a condition on L150 that it 
        # is at least 2 levels above lowest level
        #--------------------------------------------
        for c in np.arange( ncol ):
            L150[c] = np.minimum( [ L150[c] ] , [nz-3]   )[0]
        for c in np.arange( ncol ):
            te150[c] = te[L150[c], c ]
            pm150[c] = pmid[L150[c], c ]

                
    if ( Gridkey == 'tzc' ):
        nt, nz, ncol = np.shape( pmid )
        # Set up necessary arrays
        zagl  = np.zeros( (nt, search_up_L+1, ncol) )
        dzagl = np.zeros( (nt, search_up_L, ncol) )
        zaglm = np.zeros( (nt, search_up_L, ncol) )
        k150x = np.zeros( (nt, search_up_L, ncol) , dtype=np.int )
        te150 = np.zeros( (nt, ncol) )
        pm150 = np.zeros( (nt, ncol) )
        for i in np.arange( nt ):
            for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
                L = nz - 1 + (Lx - (search_up_L-1) )
                dzagl[i,Lx,:] = ( Rgas *te[i,L,:]/(grav*pmid[i,L,:]) ) * delp[i,L,:]
            for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
                zagl[i,Lx,:] = zagl[i,Lx+1,:] + dzagl[i,Lx,:] 
            for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
                zaglm[i,Lx,:] = 0.5*(zagl[i,Lx+1,:] + zagl[i,Lx,:] )

        # Insane indexing gymnastics for 'efficiency'.
        # First line is straightforward
        k150x = np.where( zaglm < 150. , 1, 0 )
        k150  = search_up_L - np.sum( k150x , axis=1 ) -1
        L150  = nz - 1 + (k150 - (search_up_L-1) )

        # Enforce a condition on L150 that it 
        # is at least 2 levels above lowest level
        #--------------------------------------------
        for i in np.arange( nt ):
            for c in np.arange( ncol ):
                L150[i,c] = np.minimum( [ L150[i,c] ] , [nz-3]   )[0]
        for i in np.arange( nt ):
            for c in np.arange( ncol ):
                te150[i,c] = te[i, L150[i,c], c ]
                pm150[i,c] = pmid[i, L150[i,c], c ]
        
    if ( Gridkey == 'tzyx' ):
        nt, nz, ny,nx = np.shape( pmid )
        # Set up necessary arrays
        zagl  = np.zeros( (nt, search_up_L+1, ny,nx) )
        dzagl = np.zeros( (nt, search_up_L, ny,nx) )
        zaglm = np.zeros( (nt, search_up_L, ny,nx) )
        k150x = np.zeros( (nt, search_up_L, ny,nx) , dtype=np.int )
        te150 = np.zeros( (nt, ny,nx) )
        pm150 = np.zeros( (nt, ny,nx) )
        for i in np.arange( nt ):
            for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
                L = nz - 1 + (Lx - (search_up_L-1) )
                dzagl[i,Lx,:,:] = ( Rgas *te[i,L,:,:]/(grav*pmid[i,L,:,:]) ) * delp[i,L,:,:]
            for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
                zagl[i,Lx,:,:] = zagl[i,Lx+1,:,:] + dzagl[i,Lx,:,:] 
            for Lx in np.arange( start=search_up_L-1, stop=-1, step=-1 ):
                zaglm[i,Lx,:,:] = 0.5*(zagl[i,Lx+1,:,:] + zagl[i,Lx,:,:] )

        # Insane indexing gymnastics for 'efficiency'.
        # First line is straightforward
        k150x = np.where( zaglm < 150. , 1, 0 )
        k150  = search_up_L - np.sum( k150x , axis=1 ) -1
        L150  = nz - 1 + (k150 - (search_up_L-1) )

        # Enforce a condition on L150 that it 
        # is at least 2 levels above lowest level
        #--------------------------------------------
        for i in np.arange( nt ):
            for y in np.arange( ny ):
                for x in np.arange( nx ):
                    L150[i,y,x] = np.minimum( [ L150[i,y,x] ] , [nz-3]   )[0]
        for i in np.arange( nt ):
            for y in np.arange( ny ):
                for x in np.arange( nx ):
                    te150[i,y,x] = te[i, L150[i,y,x], y,x ]
                    pm150[i,y,x] = pmid[i, L150[i,y,x], y,x ]

    

    return te150,pm150,L150

def Pressure_TandP150  ( am, bm, ai, bi, ps ,te , p_00=100_000., Gridkey='tzc' ):

    pmid, pint, delp = Pressure ( am=am ,
                                      bm=bm ,
                                      ai=ai ,
                                      bi=bi ,
                                      ps=ps ,
                                      p_00=p_00, 
                                      Gridkey=Gridkey )
    te150 , pm150 , L150  = TandP150 ( te=te ,
                                   pmid = pmid  , 
                                   delp = delp  , 
                                   search_up_L=10 , 
                                   Gridkey=Gridkey )    

    return te150,pm150,L150
