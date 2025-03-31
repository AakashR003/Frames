
try:
    from config import config
except:
    from .config import config

class NeumanBC():
    
    def __init__(self,**kwargs):
        self.type = kwargs.get("type", None)
        self.Magnitude = kwargs.get("Magnitude", None)
        self.Distance1 = kwargs.get("Distance1", None)
        self.Distance2 = kwargs.get("Distance2", None)
        self.AssignedTo = kwargs.get("AssignedTo", None)
        self.Members = kwargs.get("Members", None)
        
        self.MemberNo = int(self.AssignedTo.split()[1])-1
    
    def EquivalentLoad(self):
        
        FEDivision = config.get_FEDivision()
        self.frml=[] # Free moment Distribution(Simply supported) along beam 
        tarea=0
        tyda=0
        mp=0
        if self.type == "PL" :
            va=-self.Magnitude*(self.Members[self.MemberNo].length()-self.Distance1)/(self.Members[self.MemberNo].length())
            vb=-self.Magnitude*self.Distance1/(self.Members[self.MemberNo].length())
            while mp<=self.Members[self.MemberNo].length():
                if mp > self.Distance1 :
                    mppl=(self.Magnitude)*(mp-self.Distance1)
                else:
                    mppl=0
                m=va*mp+mppl
                area=m*(self.Members[self.MemberNo].length()/FEDivision)
                yda=area*mp
                tarea=area+tarea
                tyda=yda+tyda
                self.frml.append(m)
                mp=mp+(self.Members[self.MemberNo].length()/FEDivision)
            if(tarea==0):
                centroid=0
            else:
                centroid=tyda/tarea
            
        elif self.type == "UDL" :
            self.Range = abs(self.Distance2 - self.Distance1)
            va=-self.Magnitude*self.Range*(self.Members[self.MemberNo].length()-self.Distance1-self.Range*0.5)/(self.Members[self.MemberNo].length())
            vb=-self.Magnitude*self.Range*(self.Distance1+self.Range*0.5)/(self.Members[self.MemberNo].length())
            while mp<=self.Members[self.MemberNo].length():
                if(mp>self.Distance1 and mp<=(self.Distance1+self.Range)):
                    mpu=self.Magnitude*0.5*(mp-self.Distance1)**2
                elif(mp>(self.Distance1+self.Range) and mp<self.Members[self.MemberNo].length()):
                    mpu=self.Magnitude*self.Range*(self.Range*0.5+(mp-(self.Distance1+self.Range)))
                else:
                    mpu=0
                m=va*mp+mpu
                area=m*(self.Members[self.MemberNo].length()/FEDivision)
                yda=area*mp
                tarea=area+tarea
                tyda=yda+tyda
                self.frml.append(m)
                mp=mp+(self.Members[self.MemberNo].length()/FEDivision)
            if(tarea==0):
                centroid=0
            else:
                centroid=tyda/tarea
        
        else:
            raise ValueError(f"Unsupported load type: '{self.type}'")
        
        self.mfab=((2*(tarea*(self.Members[self.MemberNo].length()-centroid)*6/self.Members[self.MemberNo].length()/self.Members[self.MemberNo].length())-(tarea*centroid*6/self.Members[self.MemberNo].length()/self.Members[self.MemberNo].length()))/3)
        self.mfba=-(2*(tarea*centroid*6/self.Members[self.MemberNo].length()/self.Members[self.MemberNo].length())-(tarea*(self.Members[self.MemberNo].length()-centroid)*6/self.Members[self.MemberNo].length()/self.Members[self.MemberNo].length()))/3

        
        if(self.Members[self.MemberNo].alpha()>=0):
            self.V_b=(-self.mfab-self.mfba+vb*self.Members[self.MemberNo].length())/self.Members[self.MemberNo].length()
            self.V_a=(self.mfab+self.mfba+va*self.Members[self.MemberNo].length())/self.Members[self.MemberNo].length()
        else:
            self.V_b=-(-self.mfab-self.mfba+vb*self.Members[self.MemberNo].length())/self.Members[self.MemberNo].length()
            self.V_a=-(self.mfab+self.mfba+va*self.Members[self.MemberNo].length())/self.Members[self.MemberNo].length()
        
        
        return {"Ha":(0,self.Members[self.MemberNo].DoFNumber()[0]),
                "Va":(self.V_a,self.Members[self.MemberNo].DoFNumber()[1]),
                "Ma":(self.mfab,self.Members[self.MemberNo].DoFNumber()[2]),
                "Hb":(0,self.Members[self.MemberNo].DoFNumber()[3]),
                "Vb":(self.V_b,self.Members[self.MemberNo].DoFNumber()[4]),
                "Mb":(self.mfba,self.Members[self.MemberNo].DoFNumber()[5]),
                "FreeMoment": self.frml}#[self.V_a,self.mfab,self.V_b,self.mfba]
