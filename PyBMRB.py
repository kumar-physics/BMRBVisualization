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
    def plotn15_hsqc(self, entryids,colorby):
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
        return isOk, msg, inData

    @staticmethod
    def GetSTARInfo(starData, dflag):
        sf_list = []
        lp_list = []
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
        return sf_list, lp_list

    def getCSdata(self, starData):
        cs_loop = starData.get_loops_by_category('_Atom_chem_shift')
        out = []
        for cs in cs_loop:
            tout = []
            rid_list = sorted([int(i) for i in list(set(cs.get_tag(['Comp_index_ID'])))])
            cs_dat = cs.get_tag(['Comp_index_ID', 'Comp_ID', 'Atom_ID', 'Atom_type', 'Val'])
            cs_dict = {}
            for c in cs_dat:
                if int(c[0]) not in cs_dict.keys(): cs_dict[int(c[0])] = []
                cs_dict[int(c[0])].append(c)
            for rid in rid_list:
                cs_out = get_pair('N', 'H', rid, cs_dict[rid])
                if cs_out[2] is not None and cs_out[3] is not None:
                    tout.append(cs_out)
            out.append(tout)
        return out

    def plot(self, csdata):
        x = []
        y = []
        l = {}
        for i in csdata:
            if i[1] not in l.keys(): l[i[1]] = [[], [], []]

            l[i[1]][0].append(i[2])
            l[i[1]][1].append(i[3])
            l[i[1]][2].append('{}-{}'.format(i[0], i[1]))
        # print l
        data = []
        for k in l.keys():
            data.append(
                plotly.graph_objs.Scatter(
                    x=l[k][1],
                    y=l[k][0],
                    text=l[k][2],
                    mode='markers',
                    name=k)
            )
        #         t1= plotly.graph_objs.Scatter(x=x,y=y,mode = 'markers',text=l)
        #         data = [t1,t1]
        layout = plotly.graph_objs.Layout(
            xaxis=dict(autorange='reversed',
                       title='H (ppm)'),
            yaxis=dict(autorange='reversed',
                       title='N (ppm)'))
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        plotly.offline.plot(fig, filename='hsqc.html', auto_open=False)

    def histogram(self, atom):
        link = 'http://webapi.bmrb.wisc.edu/v2/search/chemical_shifts?atom_id={}'.format(atom)
        x = requests.get(link, headers={"Application": "BMRB-Plotly"})
        dat = x.json()['data']
        col = x.json()['columns']
        rid = col.index('Atom_chem_shift.Comp_ID')
        csid = col.index('Atom_chem_shift.Val')
        csdict = {}
        for i in dat:

            if i[rid] not in csdict.keys(): csdict[i[rid]] = []
            csdict[i[rid]].append(i[csid])

        hist_data = []
        group_lables = []
        # print csdict
        for k in csdict.keys():
            if len(csdict[k]) > 500:
                hist_data.append(csdict[k])
                group_lables.append(k)
                #print k, len(csdict[k])
        fig = plotly.tools.FigureFactory.create_distplot(hist_data, group_lables, bin_size=0.5)
        plotly.offline.plot(fig, filename='test2')


def readBMRBentry(entryid):
    isOk = False

    try:
        inData = pynmrstar.Entry.from_database(entryid)
        isOk = True
    except ValueError:
        print ("problem")

    return (isOk, inData)


def get_pair(atom1, atom2, resid, cslist):
    v1 = None
    v2 = None
    for x in cslist:
        if atom1 == x[2]:
            v1 = x[4]
        if atom2 == x[2]:
            v2 = x[4]

    return [resid, x[1], v1, v2]


if __name__ == "__main__":
    # id = sys.argv[1]
    # p = BMRBVisualization()
    # # p.histogram('HD2')
    # x = readBMRBentry(id)
    # z = p.getCSdata(x[1])[0]
    # p.plot(z)
    p = PyBMRB()
