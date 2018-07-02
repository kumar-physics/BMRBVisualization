from __future__ import print_function
import pynmrstar
import plotly
from __builtin__ import staticmethod
import requests
import sys



class PyBMRB(object):

    def __init__(self):
        print ("PyNMRSTAR version : {}".format(pynmrstar.__version__))
        #self.plotn15_hsqc([17074,17076,17077],'entry')

    def getEntry(self,entryid):
        if type(entryid) is list:
            outdata=[]
            for eid in entryid:
                indata = pynmrstar.Entry.from_database(eid)
                cs_data=indata.get_tags(['_Atom_chem_shift.Comp_index_ID', '_Atom_chem_shift.Comp_ID', '_Atom_chem_shift.Atom_ID',
                    '_Atom_chem_shift.Atom_type', '_Atom_chem_shift.Assigned_chem_shift_list_ID','_Atom_chem_shift.Val'])

                eids = [eid for i in range(len(cs_data['_Atom_chem_shift.Comp_index_ID']))]
                eid_cs_data=[eids,cs_data['_Atom_chem_shift.Comp_index_ID'],
                             cs_data['_Atom_chem_shift.Comp_ID'],
                             cs_data['_Atom_chem_shift.Atom_ID'],
                             cs_data['_Atom_chem_shift.Atom_type'],
                             cs_data['_Atom_chem_shift.Assigned_chem_shift_list_ID'],
                             cs_data['_Atom_chem_shift.Val']]
                if len(outdata):
                    for i in range(len(eid_cs_data)):
                        outdata[i] = outdata[i] + eid_cs_data[i]
                else:
                    outdata = eid_cs_data
        else:
            indata = pynmrstar.Entry.from_database(entryid)
            cs_data = indata.get_tags(
                ['_Atom_chem_shift.Comp_index_ID', '_Atom_chem_shift.Comp_ID', '_Atom_chem_shift.Atom_ID',
                 '_Atom_chem_shift.Atom_type', '_Atom_chem_shift.Assigned_chem_shift_list_ID', '_Atom_chem_shift.Val'])
            eids = [entryid for i in range(len(cs_data['_Atom_chem_shift.Comp_index_ID']))]
            outdata = [eids,cs_data['_Atom_chem_shift.Comp_index_ID'],
                           cs_data['_Atom_chem_shift.Comp_ID'],
                           cs_data['_Atom_chem_shift.Atom_ID'],
                           cs_data['_Atom_chem_shift.Atom_type'],
                           cs_data['_Atom_chem_shift.Assigned_chem_shift_list_ID'],
                       cs_data['_Atom_chem_shift.Val']]

        return outdata

    def n15_hsqc(self,csdata):
        sidechainres = ['ARG','GLN','ASN','HIS','TRP','LYS']
        sidechains = {
            'ARG-HH11':['HH11''NH1'],
            'ARG-HH12':['HH12', 'NH1'],
            'ARG-HH21':['HH21','NH2'],
            'ARG-HH22':['HH22', 'NH2'],
            'ARG-HE':['HE','NE'],
            'GLN-HE21':['HE21','NE2'],
            'GLN-HE22':['HE22', 'NE2'],
            'ASN-HD21':['HD21','ND2'],
            'ASN-HD22':['HD22', 'ND2'],
            'HIS-HD1':['HD1''ND1'],
            'HIS-HE2':['HE2','NE2'],
            'TRP-HE1':['HE1','NE1'],
            'LYS-HZ':['HZ','NZ'],
            'LYS-HZ1': ['HZ1',  'NZ'],
            'LYS-HZ2': ['HZ2',  'NZ'],
            'LYS-HZ3': [ 'HZ3', 'NZ']
        }
        outdata = [[],[],[],[],[]]
        for i in range(len(csdata[0])):
            atomid = '{}-{}-{}-{}'.format(csdata[0][i],csdata[1][i],csdata[2][i],csdata[5][i])
            if csdata[3][i]=="H":
                if atomid not in outdata[0]:
                    outdata[0].append(atomid)
                    outdata[1].append(csdata[6][i])
                    outdata[2].append(None)
                    outdata[3].append(csdata[3][i])
                    outdata[4].append(None)

                else:
                    outdata[1][outdata[0].index(atomid)]=csdata[6][i]
                    outdata[3][outdata[0].index(atomid)]=csdata[3][i]
            if csdata[3][i]=="N":
                if atomid not in outdata[0]:
                    outdata[0].append(atomid)
                    outdata[2].append(csdata[6][i])
                    outdata[1].append(None)
                    outdata[4].append(csdata[3][i])
                    outdata[3].append(None)
                else:
                    outdata[2][outdata[0].index(atomid)] = csdata[6][i]
                    outdata[4][outdata[0].index(atomid)]=(csdata[3][i])
            if csdata[2][i] in sidechainres:
                for k in sidechains.keys():
                    if k.split("-")[0]== csdata[2][i]:
                        atomid = '{}-{}-{}-{}'.format(csdata[0][i],csdata[1][i],k,csdata[5][i])
                        if csdata[3][i] in sidechains[k] and csdata[4][i]== "H":
                            if atomid not in outdata[0]:
                                outdata[0].append(atomid)
                                outdata[1].append(csdata[6][i])
                                outdata[2].append(None)
                                outdata[3].append(csdata[3][i])
                                outdata[4].append(None)

                            else:
                                outdata[1][outdata[0].index(atomid)] = csdata[6][i]
                                outdata[3][outdata[0].index(atomid)] = csdata[3][i]
                        if csdata[3][i] in sidechains[k] and csdata[4][i] == "N":
                            if atomid not in outdata[0]:
                                outdata[0].append(atomid)
                                outdata[2].append(csdata[6][i])
                                outdata[1].append(None)
                                outdata[4].append(csdata[3][i])
                                outdata[3].append(None)
                            else:
                                outdata[2][outdata[0].index(atomid)] = csdata[6][i]
                                outdata[4][outdata[0].index(atomid)] = (csdata[3][i])

        return outdata
    def plotn15_hsqc(self, entryids,colorby= 'res'):
        csdata = self.getEntry(entryids)
        hsqcdata = self.n15_hsqc(csdata)
        if colorby == 'entry':
            id=0
        elif colorby == 'res':
            id = 2
        else:
            print("Colorby error")

        groups = set([k.split("-")[id] for k in hsqcdata[0]])
        data_sets = {}
        for gid in groups:
            data_sets[gid]=[[],[],[]]
            for i in range(len(hsqcdata[0])):
                if hsqcdata[0][i].split("-")[id] == gid:
                    data_sets[gid][0].append(hsqcdata[1][i])
                    data_sets[gid][1].append(hsqcdata[2][i])
                    data_sets[gid][2].append(hsqcdata[0][i])
        data = []
        for k in data_sets.keys():
            data.append(plotly.graph_objs.Scatter(x=data_sets[k][0],
                                          y=data_sets[k][1],
                                                  text = data_sets[k][2],
                                          mode='markers',
                                          name=k)
                                          )
        layout = plotly.graph_objs.Layout(
            xaxis=dict(autorange='reversed',
                       title='H (ppm)'),
            yaxis=dict(autorange='reversed',
                       title='N (ppm)'),
            showlegend=True,
            hovermode='closest')
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        plotly.offline.plot(fig, filename='hsqc.html', auto_open=True)






if __name__ == "__main__":
    p = PyBMRB()
