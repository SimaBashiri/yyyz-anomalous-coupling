from HiggsAnalysis.CombinedLimit.PhysicsModel import *



class YYYZ(PhysicsModelBase):
    def __init__(self):
        self.options = {}
        self.verbose = False


    def setPhysicsOptions(self, physOptions):
        for po in physOptions:
            if po.startswith("POI="):
                self.poi = po.split("=")[1]

    def doParametersOfInterest(self):
        """Create POI and other parameters, and define the POI set."""
        self.modelBuilder.doVar("zta[1e-1,0.0,1]")
        self.modelBuilder.doVar("zta~[1e-2,0.0,1]") 
        self.modelBuilder.doVar("r[1,0.0,10]") 
        self.modelBuilder.factory_('expr::term1("((3*@0*@0)+(3*@1*@1)-(2*@0*@1))*@2/0.0283", zta,zta~,r)') 
        self.modelBuilder.doSet("POI", "zta,zta~,r")



    def getYieldScale(self, bin, process):

        if process == "signal":
            return 'term1'

        else:
            return 1 


yyyz = YYYZ()

