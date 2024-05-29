class VariableContainer:
    def __init__(self):
        self.MyDst=None
        self.MyDstVgrid=None
        self.MySrc=None

        self.regrd=None
        self.srcf=None
        self.dstf=None

        self.phis_CAM=None
        self.phis_ERA=None

        self.amid_CAM=None
        self.bmid_CAM=None
        self.aint_CAM=None
        self.bint_CAM=None

        self.lon_CAM=None
        self.lat_CAM=None
        self.area_CAM=None
        self.area_ERA=None

        # Grid keys for remapping
        self.srcHkey=None
        self.dstHkey=None
        self.srcTHkey=None
        self.dstTHkey=None
        self.srcZHkey=None
        self.dstZHkey=None
        self.srcTZHkey=None
        self.dstTZHkey=None

        self.doWilliamsonOlson=None
        self.RegridMethod=None

        self.p_00_ERA=None
        self.p_00_CAM=None

        self.dstHkey = None
        self.dst_type=None
        self.dst_scrip = None 
        self.dst_TopoFile = None 
        self.dstVgridFile=None

        self.srcHkey = None
        self.src_type= None
        self.src_scrip = None 
        self.src_TopoFile = None


        self.pdTime_ERA = None
        self.ps_ERA = None
        self.te_ERA = None
        self.q_ERA = None
        self.u_ERA = None
        self.v_ERA = None
        self.w_ERA = None
        self.amid_ERA  = None
        self.bmid_ERA   = None
        self.aint_ERA   = None
        self.bint_ERA = None

        # For diagnostic puroposes
        self.lon_ERA = None
        self.lat_ERA = None


        self.pmid_ERA = None 
        self.pint_ERA = None 
        self.delp_ERA = None 
        self.pmid_CAM_zERA  = None 
        self.pint_CAM_zERA  = None 
        self.delp_CAM_zERA  = None 
        self.pmid_CAM = None 
        self.pint_CAM = None 
        self.delp_CAM  = None 

        self.te_150 = None
        self.pmid_150 = None
        self.L150 = None

        # Variables created by Regridder
        #----------------------------------------------------------
        self.ps_CAM   = None 
        self.ps_ERA_xCAM = None
        self.phis_ERA_xCAM = None
        self.te_ERA_xzCAM = None
        self.q_ERA_xzCAM = None
        self.u_ERA_xzCAM = None
        self.v_ERA_xzCAM = None
        self.w_ERA_xzCAM = None

        self.te_ERA_xCAM_00 = None
        self.te_ERA_xCAM = None
        self.q_ERA_xCAM = None
        self.u_ERA_xCAM = None
        self.v_ERA_xCAM = None
        self.w_ERA_xCAM = None

        self.te_WO = None
        self.ts_extrap = None


        
Gv = VariableContainer()