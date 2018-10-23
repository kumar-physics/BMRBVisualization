from __future__ import print_function
import pynmrstar
import plotly
import sys
import csv
import numpy as np
import json

# Determine if we are running in python3
PY3 = (sys.version_info[0] == 3)

# pylint: disable=wrong-import-position,no-name-in-module
# pylint: disable=import-error,wrong-import-order
# Python version dependent loads
if PY3:
    from urllib.request import urlopen, Request
    from urllib.error import HTTPError, URLError
    from io import StringIO, BytesIO
else:
    from urllib2 import urlopen, HTTPError, URLError, Request
    from cStringIO import StringIO

    BytesIO = StringIO

_API_URL = "http://webapi.bmrb.wisc.edu/v2"
_NOTEBOOK = False

# http://webapi.bmrb.wisc.edu/v2/search/chemical_shifts?comp_id=ASP&atom_id=HD2

class Spectra(object):

    def __init__(self):
        print("PyNMRSTAR version : {}".format(pynmrstar.__version__))
        if _NOTEBOOK:
            plotly.offline.init_notebook_mode(connected=True)
        # self.plotn15_hsqc([17074,17076,17077],'entry')

    def getEntry(self, entryid):
        if type(entryid) is list:
            outdata = []
            for eid in entryid:
                indata = pynmrstar.Entry.from_database(eid)
                cs_data = indata.get_tags(
                    ['_Atom_chem_shift.Comp_index_ID', '_Atom_chem_shift.Comp_ID', '_Atom_chem_shift.Atom_ID',
                     '_Atom_chem_shift.Atom_type', '_Atom_chem_shift.Assigned_chem_shift_list_ID',
                     '_Atom_chem_shift.Val'])

                eids = [eid for i in range(len(cs_data['_Atom_chem_shift.Comp_index_ID']))]
                eid_cs_data = [eids, cs_data['_Atom_chem_shift.Comp_index_ID'],
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
            outdata = [eids, cs_data['_Atom_chem_shift.Comp_index_ID'],
                       cs_data['_Atom_chem_shift.Comp_ID'],
                       cs_data['_Atom_chem_shift.Atom_ID'],
                       cs_data['_Atom_chem_shift.Atom_type'],
                       cs_data['_Atom_chem_shift.Assigned_chem_shift_list_ID'],
                       cs_data['_Atom_chem_shift.Val']]

        return outdata

    def n15_hsqc(self, csdata):
        sidechainres = ['ARG', 'GLN', 'ASN', 'HIS', 'TRP', 'LYS']
        sidechains = {
            'ARG-HH11': ['HH11''NH1'],
            'ARG-HH12': ['HH12', 'NH1'],
            'ARG-HH21': ['HH21', 'NH2'],
            'ARG-HH22': ['HH22', 'NH2'],
            'ARG-HE': ['HE', 'NE'],
            'GLN-HE21': ['HE21', 'NE2'],
            'GLN-HE22': ['HE22', 'NE2'],
            'ASN-HD21': ['HD21', 'ND2'],
            'ASN-HD22': ['HD22', 'ND2'],
            'HIS-HD1': ['HD1''ND1'],
            'HIS-HE2': ['HE2', 'NE2'],
            'TRP-HE1': ['HE1', 'NE1'],
            'LYS-HZ': ['HZ', 'NZ'],
            'LYS-HZ1': ['HZ1', 'NZ'],
            'LYS-HZ2': ['HZ2', 'NZ'],
            'LYS-HZ3': ['HZ3', 'NZ']
        }
        outdata = [[], [], [], [], []]
        for i in range(len(csdata[0])):
            atomid = '{}-{}-{}-{}'.format(csdata[0][i], csdata[1][i], csdata[2][i], csdata[5][i])
            if csdata[3][i] == "H":
                if atomid not in outdata[0]:
                    outdata[0].append(atomid)
                    outdata[1].append(csdata[6][i])
                    outdata[2].append(None)
                    outdata[3].append(csdata[3][i])
                    outdata[4].append(None)

                else:
                    outdata[1][outdata[0].index(atomid)] = csdata[6][i]
                    outdata[3][outdata[0].index(atomid)] = csdata[3][i]
            if csdata[3][i] == "N":
                if atomid not in outdata[0]:
                    outdata[0].append(atomid)
                    outdata[2].append(csdata[6][i])
                    outdata[1].append(None)
                    outdata[4].append(csdata[3][i])
                    outdata[3].append(None)
                else:
                    outdata[2][outdata[0].index(atomid)] = csdata[6][i]
                    outdata[4][outdata[0].index(atomid)] = (csdata[3][i])
            if csdata[2][i] in sidechainres:
                for k in sidechains.keys():
                    if k.split("-")[0] == csdata[2][i]:
                        atomid = '{}-{}-{}-{}'.format(csdata[0][i], csdata[1][i], k, csdata[5][i])
                        if csdata[3][i] in sidechains[k] and csdata[4][i] == "H":
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

    def plotn15_hsqc(self, entryids, colorby='res', groupbyres=False):
        if type(entryids) is list:
            outfilename = 'n15hsqc.html'
            title = 'Simulated N15-HSQC peak positions'
        else:
            outfilename = '{}.html'.format(entryids)
            title = 'Simulated N15-HSQC peak positions of BMRB entry {}'.format(entryids)
        csdata = self.getEntry(entryids)
        hsqcdata = self.n15_hsqc(csdata)
        if colorby == 'entry':
            id = 0
        elif colorby == 'res':
            id = 2
        else:
            print("Colorby error")

        groups = set([k.split("-")[id] for k in hsqcdata[0]])
        data_sets = {}
        for gid in groups:
            data_sets[gid] = [[], [], []]
            for i in range(len(hsqcdata[0])):
                if hsqcdata[0][i].split("-")[id] == gid:
                    data_sets[gid][0].append(hsqcdata[1][i])
                    data_sets[gid][1].append(hsqcdata[2][i])
                    data_sets[gid][2].append(hsqcdata[0][i])

        if groupbyres:
            id = 1
            groups2 = set(["-".join(k.split("-")[1:4]) for k in hsqcdata[0]])
            data_sets2 = {}
            for gid in groups2:
                data_sets2[gid] = [[], [], []]
                for i in range(len(hsqcdata[0])):
                    if "-".join(hsqcdata[0][i].split("-")[1:4]) == gid:
                        data_sets2[gid][0].append(hsqcdata[1][i])
                        data_sets2[gid][1].append(hsqcdata[2][i])
                        data_sets2[gid][2].append(hsqcdata[0][i])

        data = []
        for k in data_sets.keys():
            data.append(plotly.graph_objs.Scatter(x=data_sets[k][0],
                                                  y=data_sets[k][1],
                                                  text=data_sets[k][2],
                                                  mode='markers',
                                                  name=k)
                        )
        if groupbyres:
            for k in data_sets2.keys():
                data.append(plotly.graph_objs.Scatter(x=data_sets2[k][0],
                                                      y=data_sets2[k][1],
                                                      text=data_sets2[k][2],
                                                      mode='lines',
                                                      name=k,
                                                      showlegend=False)
                            )

        layout = plotly.graph_objs.Layout(
            xaxis=dict(autorange='reversed',
                       title='H (ppm)'),
            yaxis=dict(autorange='reversed',
                       title='N (ppm)'),
            showlegend=True,
            hovermode='closest',
            title=title)
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        if _NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=outfilename, auto_open=True)


class Histogram(object):

    def __init__(self):
        self.data_dir = '/home/kumaran/bmrbvis'
        if _NOTEBOOK:
            plotly.offline.init_notebook_mode(connected=True)

    def get_histogram_api(self, residue, atom, filtered=True, sd_limit=10, normalized=False):
        url = Request(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue, atom))
        url.add_header('Application', 'BMRBiViz')
        r = urlopen(url)
        d = json.loads(r.read())
        x = [i[d['columns'].index('Atom_chem_shift.Val')] for i in d['data']]
        if filtered:
            mean = np.mean(x)
            sd = np.std(x)
            lb = mean - (sd_limit * sd)
            ub = mean + (sd_limit * sd)
            x = [i for i in x if lb < i and i  < ub]
        if normalized:
            data = plotly.graph_objs.Histogram(x=x, name="{}-{}".format(residue, atom), histnorm='probability',opacity=0.75)
        else:
            data = plotly.graph_objs.Histogram(x=x, name="{}-{}".format(residue, atom),opacity=0.75)
        return data

    def get_conditional_histogram_api(self, residue, atom, atomlist, cslist, filtered=True, sd_limit=10, normalized=False):
        url =Request( _API_URL + "/search/chemical_shifts?comp_id={}".format(residue))
        url.add_header('Application','BMRBiViz')
        r = urlopen(url)
        d1 = json.loads(r.read())
        d={}
        entry_id_index = d1['columns'].index('Atom_chem_shift.Entry_ID')
        seq_id_index = d1['columns'].index('Atom_chem_shift.Comp_index_ID')
        res_id_index = d1['columns'].index('Atom_chem_shift.Comp_ID')
        atom_id_index= d1['columns'].index('Atom_chem_shift.Atom_ID')
        cs_id_index = d1['columns'].index('Atom_chem_shift.Val')
        for i in d1['data']:
            entry_id = i[entry_id_index]
            seq_id = i[seq_id_index]
            res_id = i[res_id_index]
            atom_id = i[atom_id_index]
            d['{}-{}-{}-{}'.format(entry_id,seq_id,res_id,atom_id)] = i[cs_id_index]
        filter_list = []
        for k in d.keys():
            for i in range(len(atomlist)):
                if 'H' in atomlist[i]:
                    epsilon = 0.1
                else:
                    epsilon = 0.5
                if k.split("-")[-1]==atomlist[i] and (cslist[i] + epsilon < d[k] or d[k] < cslist[i] - epsilon):
                    filter_list.append('{}-{}'.format(k.split("-")[0],k.split("-")[1]))
        x=[]
        filter_list = list(set(filter_list))
        # for i in filter_list:
        #     try:
        #         k='{}-{}-{}-{}'.format(i.split("-")[0],i.split("-")[1],residue,atom)
        #         x.append(d[k])
        #     except KeyError:
        #         pass
        for k in d.keys():
            atm_id = '{}-{}'.format(k.split("-")[0],k.split("-")[1])
            if atm_id not in filter_list and k.split("-")[2] == residue and k.split("-")[3]== atom:
                #ent_id = '{}-{}-{}-{}'.format(k.split("-")[0],k.split("-")[1],residue,atom)
                x.append(d[k])


        if filtered:
            mean = np.mean(x)
            sd = np.std(x)
            lb = mean - (sd_limit * sd)
            ub = mean + (sd_limit * sd)
            x = [i for i in x if i > lb and i < ub]
        filter_values = ''
        for i in range(len(atomlist)):
            if i==len(atomlist)-1:
                filter_values += '{}:{}'.format(atomlist[i], cslist[i])
            else:
                filter_values += '{}:{},'.format(atomlist[i],cslist[i])
        if normalized:
            data = plotly.graph_objs.Histogram(x=x, name="{}-{}({})".format(residue, atom, filter_values), histnorm='probability', opacity=0.75)
        else:
            data = plotly.graph_objs.Histogram(x=x, name="{}-{}({})".format(residue, atom,filter_values, opacity=0.75))
        return data





    def get_histogram2d_api(self, residue1, atom1, residue2, atom2, filtered=True, sd_limit=10, normalized=False):
        url1 = Request(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue1, atom1))
        url2 = Request(_API_URL + "/search/chemical_shifts?comp_id={}&atom_id={}".format(residue2, atom2))
        url1.add_header('Application', 'BMRBiViz')
        url2.add_header('Application', 'BMRBiViz')
        r1 = urlopen(url1)
        r2 = urlopen(url2)
        d1 = json.loads(r1.read())
        d = {}
        for i in d1['data']:
            entry_id = i[d1['columns'].index('Atom_chem_shift.Entry_ID')]
            seq_id = i[d1['columns'].index('Atom_chem_shift.Comp_index_ID')]
            d["{}-{}".format(entry_id, seq_id)] = i[d1['columns'].index('Atom_chem_shift.Val')]
        # x = [i[d1['columns'].index('Atom_chem_shift.Val')] for i in d1['data']]
        d2 = json.loads(r2.read())
        x = []
        y = []
        for i in d2['data']:
            entry_id = i[d2['columns'].index('Atom_chem_shift.Entry_ID')]
            seq_id = i[d2['columns'].index('Atom_chem_shift.Comp_index_ID')]
            try:
                k = "{}-{}".format(entry_id, seq_id)
                x.append(d[k])
                y.append(i[d2['columns'].index('Atom_chem_shift.Val')])
            except KeyError:
                pass

        # y = [i[d2['columns'].index('Atom_chem_shift.Val')] for i in d2['data']]
        if filtered:
            meanx = np.mean(x)
            meany = np.mean(y)
            sdx = np.std(x)
            sdy = np.std(y)
            lbx = meanx - (sd_limit * sdx)
            lby = meany - (sd_limit * sdy)
            ubx = meanx + (sd_limit * sdx)
            uby = meany + (sd_limit * sdy)
            x = [i for i in x if i > lbx and i < ubx]
            y = [i for i in y if i > lby and i < uby]

            if 'H' in atom1:
                binsizex = 0.01
            else:
                binsizex = 0.5

            if 'H' in atom2:
                binsizey = 0.01
            else:
                binsizey = 0.5

            nbinsx = round((max(x)-min(x))/binsizex)
            nbinsy = round((max(y) - min(y)) / binsizey)
            #print (nbinsx,nbinsy,binsizex,binsizey)

        if normalized:
            data = [plotly.graph_objs.Histogram2dContour(x=x, y=y, histnorm='probability', colorscale='Jet'),
                    plotly.graph_objs.Histogram(
                        y=y,
                        xaxis='x2',
                        name="{}-{}".format(residue2, atom2),
                        histnorm='probability'
                    ),
                    plotly.graph_objs.Histogram(
                        x=x,
                        yaxis='y2',
                        name="{}-{}".format(residue1, atom1),
                        histnorm='probability'
                    )
                    ]
        else:
            data = [plotly.graph_objs.Histogram2dContour(x=x, y=y, colorscale='Jet'),
                    plotly.graph_objs.Histogram(
                        y=y,
                        xaxis='x2',
                        nbinsy = nbinsx,
                        name="{}-{}".format(residue2, atom2)
                    ),
                    plotly.graph_objs.Histogram(
                        x=x,
                        nbinsx = nbinsy,
                        yaxis='y2',
                        name="{}-{}".format(residue1, atom1)
                    )
                    ]

        return data

    def get_histogram(self, residue, atom, filtered=True, sd_limit=10, normalized=False):
        data_file = '{}/{}_{}_sel.txt'.format(self.data_dir, residue, atom)
        with open(data_file) as csv_file:
            csv_reader = csv.reader(csv_file, delimiter=',')
            x = []
            for row in csv_reader:
                x.append(float(row[7]))
        if filtered:
            mean = np.mean(x)
            sd = np.std(x)
            lb = mean - (sd_limit * sd)
            ub = mean + (sd_limit * sd)
            x = [i for i in x if i > lb and i < ub]
        if normalized:
            data = plotly.graph_objs.Histogram(x=x, name="{}-{}".format(residue, atom), histnorm='probability')
        else:
            data = plotly.graph_objs.Histogram(x=x, name="{}-{}".format(residue, atom))
        return data

    def get_histogram2d(self, residue1, atom1, residue2, atom2, filtered=True, sd_limit=10, normalized=False):
        data_file1 = '{}/{}_{}_sel.txt'.format(self.data_dir, residue1, atom1)
        data_file2 = '{}/{}_{}_sel.txt'.format(self.data_dir, residue2, atom2)
        x = []
        y = []
        with open(data_file1) as csv_file1:
            csv_reader1 = csv.reader(csv_file1, delimiter=',')
            d1 = {}
            for row1 in csv_reader1:
                d1["-".join([row1[0], row1[3]])] = float(row1[7])
        with open(data_file2) as csv_file2:
            csv_reader2 = csv.reader(csv_file2, delimiter=',')
            for row2 in csv_reader2:
                try:
                    x.append(d1["-".join([row2[0], row2[3]])])
                    y.append(float(row2[7]))
                except KeyError:
                    pass

        if filtered:
            meanx = np.mean(x)
            meany = np.mean(y)
            sdx = np.std(x)
            sdy = np.std(y)
            lbx = meanx - (sd_limit * sdx)
            lby = meany - (sd_limit * sdy)
            ubx = meanx + (sd_limit * sdx)
            uby = meany + (sd_limit * sdy)
            x = [i for i in x if i > lbx and i < ubx]
            y = [i for i in y if i > lby and i < uby]

        if normalized:
            data = [plotly.graph_objs.Histogram2dContour(x=x, y=y, histnorm='probability'),
                    plotly.graph_objs.Histogram(
                        y=y,
                        xaxis='x2',
                        name="{}-{}".format(residue2, atom2),
                        histnorm='probability'
                    ),
                    plotly.graph_objs.Histogram(
                        x=x,
                        yaxis='y2',
                        name="{}-{}".format(residue1, atom1),
                        histnorm='probability'
                    )
                    ]
        else:
            data = [plotly.graph_objs.Histogram2dContour(x=x, y=y),
                    plotly.graph_objs.Histogram(
                        y=y,
                        xaxis='x2',
                        name="{}-{}".format(residue2, atom2)
                    ),
                    plotly.graph_objs.Histogram(
                        x=x,
                        yaxis='y2',
                        name="{}-{}".format(residue1, atom1)
                    )
                    ]
        return data

    def single_2dhistogram(self, residue, atom1, atom2, filtered=True, sd_limit=10, normalized=False):
        layout = plotly.graph_objs.Layout(
            autosize=True,
            xaxis=dict(
                zeroline=False,
                domain=[0, 0.85],
                showgrid=True,
                title='{}-{} [ppm]'.format(residue, atom1)

            ),
            yaxis=dict(
                zeroline=False,
                domain=[0, 0.85],
                showgrid=True,
                title='{}-{} [ppm]'.format(residue, atom2)

            ),
            xaxis2=dict(
                zeroline=False,
                domain=[0.85, 1],
                showgrid=True

            ),
            yaxis2=dict(
                zeroline=False,
                domain=[0.85, 1],
                showgrid=True
            ),
            hovermode='closest',
            showlegend=False
        )
        data = self.get_histogram2d_api(residue, atom1, residue, atom2, filtered, sd_limit, normalized)
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        out_file = 'histogram2d.html'
        if _NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=out_file)

    def single_atom(self, residue, atom, filtered=True, sd_limit=10, normalized=False):
        layout = plotly.graph_objs.Layout(
            barmode='overlay',
            xaxis=dict(title='Chemical Shift [ppm]'),
            yaxis=dict(title='Count'))
        data = [self.get_histogram_api(residue, atom, filtered, sd_limit, normalized)]
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        out_file = '{}_{}.html'.format(residue, atom)
        if _NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=out_file)

    def multiple_atom(self, atom_list, filtered=True, sd_limit=10, normalized=False):
        data = []
        for atm in atom_list:
            residue = atm.split("-")[0]
            atom = atm.split("-")[1]
            data.append(self.get_histogram_api(residue, atom, filtered, sd_limit, normalized))
        layout = plotly.graph_objs.Layout(
            barmode='overlay',
            xaxis=dict(title='Chemical Shift [ppm]'),
            yaxis=dict(title='Count'))
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        out_file = 'Multiple_atom_histogram.html'
        if _NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=out_file)

    def conditional_histogram(self,residue,atom, atomlist, cslist, filtered=True, sd_limit=10, normalized=False):
        layout = plotly.graph_objs.Layout(
            barmode='overlay',
            xaxis=dict(title='Chemical Shift [ppm]'),
            yaxis=dict(title='Count'))
        data = [self.get_histogram_api(residue,atom,filtered, sd_limit, normalized),
                self.get_conditional_histogram_api(residue, atom, atomlist,cslist,filtered, sd_limit, normalized)
                ]
        fig = plotly.graph_objs.Figure(data=data, layout=layout)
        out_file = '{}_{}.html'.format(residue, atom)
        if _NOTEBOOK:
            plotly.offline.iplot(fig)
        else:
            plotly.offline.plot(fig, filename=out_file)



if __name__ == "__main__":
    p = Histogram()
    #atlist = ['ASP-HA', 'GLN-HB2']
    #p.multiple_atom(atlist, normalized=False)
    p.single_2dhistogram(sys.argv[1],sys.argv[2],sys.argv[3])
    #atmlist = ['CB']
    #cslist = [69.0]
    #p.conditional_histogram('TYR','CA',atmlist,cslist)