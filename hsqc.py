import pynmrstar
import plotly
from __builtin__ import staticmethod
import requests
import grp


class BMRBVisualization(object):
    
    def __init__(self):
        pass
    
    @staticmethod 
    def readSTARFile(inFile):
        
        isOk = False
        try:
            inData = pynmrstar.Entry.from_file(inFile)
            isOk = True
            msg = "Entry"
        except ValueError as e:
            try:
                inData = pynmrstar.Saveframe.from_file(inFile)
                isOk = True
                msg = "Saveframe"
            except ValueError as e:
                try:
                    inData = pynmrstar.Loop.from_file(inFile)
                    isOk = True
                    msg = "Loop"
                except ValueError as e:
                    inData = None
                    msg = "Invalid STAR file! (or) May contain more than one Saveframe or loop in a text file "
        return (isOk,msg,inData)
    
    @staticmethod
    def GetSTARInfo(starData,dflag):
        sf_list=[]
        lp_list=[]
        if dflag == "Entry":
            for sf in starData.frame_list:
                    sf_list.append(sf.category)
                    for lp in sf:
                        lp_list.append(lp.category)
        elif dflag == "Saveframe":
            sf_list.append(starData.category)
            for lp in starData:
                lp_list.append(lp.category)
        else:
            lp_list.append(starData.category)
        return(sf_list,lp_list)
    
    
    @staticmethod
    def readBMRBentry(entryid):
        isOk = False
        
        try:
            inData = pynmrstar.Entry.from_database(entryid)
            isOk = True
        except ValueError:
            print "problem"
        
        return(isOk,inData)
    @staticmethod
    def get_pair(atom1,atom2,resid,cslist):
        v1=None
        v2=None
        for x in cslist:
            if atom1 == x[2]: 
                v1 = x[4]
            if atom2 == x[2]:
                v2 = x[4]
            
        return [resid,x[1],v1,v2]        
   
    def getCSdata(self,starData):
        cs_loop = starData.get_loops_by_category('_Atom_chem_shift')
        out=[]
        for cs in cs_loop:
            tout=[]
            rid_list= sorted([int(i) for i in list(set(cs.get_tag(['Comp_index_ID'])))])
            cs_dat=cs.get_tag(['Comp_index_ID','Comp_ID','Atom_ID','Atom_type','Val'])
            cs_dict = {}
            for c in cs_dat:
                if int(c[0]) not in cs_dict.keys(): cs_dict[int(c[0])]=[]
                cs_dict[int(c[0])].append(c) 
            for rid in rid_list:
                cs_out=self.get_pair('N','H',rid,cs_dict[rid])
                if cs_out[2] is not None and cs_out[3] is not None:
                    tout.append(cs_out)
            out.append(tout)
        return out
                
    def plot(self,csdata):
        x=[]
        y=[]
        l={}
        for i in csdata:
            if i[1] not in l.keys(): l[i[1]]=[[],[],[]]
            
            l[i[1]][0].append(i[2])
            l[i[1]][1].append(i[3])
            l[i[1]][2].append('{}-{}'.format(i[0],i[1]))
        print l
        data = []
        for k in l.keys():
            data.append(
                plotly.graph_objs.Scatter(
                    x=l[k][1],
                    y=l[k][0],
                    text = l[k][2],
                    mode = 'markers',
                    name = k)
                )
#         t1= plotly.graph_objs.Scatter(x=x,y=y,mode = 'markers',text=l)
#         data = [t1,t1]
        layout = plotly.graph_objs.Layout(
            xaxis=dict(autorange='reversed',
                       title = 'H (ppm)'),
            yaxis=dict(autorange='reversed',
                       title = 'N (ppm)'))
        fig = plotly.graph_objs.Figure(data=data,layout = layout) 
        plotly.offline.plot(fig,filename='test')
                    
    def histogram(self,atom):
        link = 'http://webapi.bmrb.wisc.edu/v2/search/chemical_shifts?atom_id={}'.format(atom)
        x=requests.get(link,headers = {"Application":"BMRB-Plotly"})
        dat= x.json()['data']
        col = x.json()['columns']
        rid = col.index('Atom_chem_shift.Comp_ID')
        csid = col.index('Atom_chem_shift.Val')
        csdict={}
        for i in dat:
           
            if i[rid] not in csdict.keys(): csdict[i[rid]]=[]
            csdict[i[rid]].append(i[csid])
            
        hist_data = []
        group_lables = []
        #print csdict
        for k in csdict.keys():
            if len(csdict[k])>500:
                hist_data.append(csdict[k])
                group_lables.append(k)
                print k,len(csdict[k])
        fig = plotly.tools.FigureFactory.create_distplot(hist_data,group_lables,bin_size=0.5)
        plotly.offline.plot(fig,filename='test2')
    
if __name__ == "__main__":
    p = BMRBVisualization()
    p.histogram('CB')
    #x=p.readBMRBentry(15007)
    #z=p.getCSdata(x[1])[0]
    #p.plot(z)
    
    
    