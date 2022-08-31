#!/usr/bin/env python
from dflow import config, s3_config
import time
import numpy as np
from dflow import (
    ShellOPTemplate,
    InputParameter,
    OutputParameter,
    InputArtifact,
    OutputArtifact,
    Workflow,
    Step,
    upload_artifact,
    download_artifact
)


from dflow.plugins.lebesgue import LebesgueContext
config['host'] = "http://39.106.93.187:32746"
s3_config["endpoint"] = "39.106.93.187:30900"


def make_gaussian_input(sys_data, fp_params):
    coordinates = sys_data['coord_t']
    atomic_symbol = sys_data['atomic_symbol']
    nproc = fp_params['nproc']; mem = fp_params['mem']
    keywords = fp_params['keywords']; pseudo = fp_params['pseudo']
    multiplicity = fp_params['multiplicity']
    charge = fp_params.get('charge',0); chk = fp_params['chk']
    chkkeywords = '%chk={}.chk'.format(chk)
    nprockeywords = '%nproc={:d}'.format(nproc)
    memkeywords = '%mem={:d}gb'.format(mem)
    titlekeywords = 'DFLOW'
    chargekeywords = '{} {}'.format(charge, multiplicity)
    buff = [chkkeywords, nprockeywords, memkeywords, '# {}'.format(keywords), '', titlekeywords, '', chargekeywords]
    for ii, (symbol, coordinate) in enumerate(zip(atomic_symbol, coordinates)):
        buff.append("%s %f %f %f" %(symbol, *coordinate))
    buff.append('')
    if pseudo == True:
        pseudo_element = fp_params['pseudo_element']
        other_element = fp_params['other_element']
        pseudo_basis = fp_params['pseudo_basis']
        other_basis = fp_params['other_basis']
        buff.append('%s 0' %(pseudo_element))
        buff.append(pseudo_basis)
        buff.append('****')
        buff.append('%s 0' %(other_element))
        buff.append(other_basis)
        buff.append('****')
        buff.append('')
        buff.append('%s 0' %(pseudo_element))
        buff.append(pseudo_basis)
        buff.append('\n')
    return '\n'.join(buff)

def s0data(file_name):
    flag = 0; energy_t = []; coord_t = []; atomic_symbol = []
    homo_t = []; lumo_t = []; lumo = None
    with open(file_name) as fp:
        for line in fp:
            if line.startswith(" Alpha  occ. eigenvalues --"):
                homo = float(line.split()[-1])*27.211386
            elif line.startswith(" Alpha virt. eigenvalues --") and lumo == None and homo != None:
                lumo = float(line.split()[4])*27.211386
                energy_t.append(energy); coord_t.append(coord); homo_t.append(homo)
                lumo_t.append(lumo); homo = None; lumo = None
            elif line.startswith(" SCF Done"):
                energy = float(line.split()[4])*27.211386
            elif line.startswith(" Symbolic Z-matrix:"): 
                flag = 1
                atomic_symbol = []
            elif line.startswith(" Center     Atomic      Atomic             Coordinates (Angstroms)"):
                flag = 4
                coord = []
            if 1 <= flag <= 2 or 4 <= flag <= 6:
                flag += 1
            elif flag == 3:
                if line.startswith(" GradGradGradGrad"):
                    flag = 0
                else:
                    s = line.strip().split()
                    if len(s) == 4:
                        atomic_symbol.append(s[0])
            elif flag == 7:
                if line.startswith(" -------"):
                    flag = 0
                else:
                    s = line.strip().split()
                    coord.append([float(x) for x in s[3:6]])
    homo_t = homo_t[-1]; lumo_t = lumo_t[-1]
    energy_t = energy_t[-1]; coord_t = coord_t[-1]
    data = {}
    data['homo_t'] = homo_t; data['lumo_t'] = lumo_t; data['energy_t'] = energy_t
    data['coord_t'] = coord_t; data['atomic_symbol'] = atomic_symbol
    return data


def s1data(file_name):
    excitation_e_t = []; energy_t = []
    with open(file_name) as fp:
        for line in fp:
            if line.startswith(" Excited State   1:      Singlet-A"):
                excitation_e_t.append(float(line.split()[6]))
            elif line.startswith(" SCF Done"):
                energy_t.append(float(line.split()[4])*27.211386)
    excitation_e_t = excitation_e_t[-1]; energy_t = energy_t[-1]
    data = {}
    data['s1_exciation_e'] = excitation_e_t
    data['s1_energy'] = energy_t
    return data

def t1data(file_name):
    energy_t = []
    with open(file_name) as fp:
        for line in fp:
            if line.startswith(" SCF Done"):
                energy_t.append(float(line.split()[4])*27.211386)
    energy_t = energy_t[-1]
    data = {}
    data['t1_energy'] = energy_t
    return data

if __name__ == "__main__":
    
    #data = s0data('s0.log')
    #fp_s1 = {'nproc':16,'mem':4,'keywords':"opt b3lyp/gen pseudo=read TD",'pseudo':True,\
    #        'pseudo_element':'Ir','other_element':'C H O N','pseudo_basis':'Lanl2DZ',\
    #        'other_basis':'6-31G**','multiplicity': 1, 'charge': 0, 'chk':'s1'}
    #ret = make_gaussian_input(data,fp_s1)
    #with open('s1.com','w') as fp:
    #    fp.write(ret)
    #fp_st = {'nproc':16,'mem':4,'keywords':"opt b3lyp/gen pseudo=read",'pseudo':True,\
    #        'pseudo_element':'Ir','other_element':'C H O N','pseudo_basis':'Lanl2DZ',\
    #        'other_basis':'6-31G**','multiplicity': 3, 'charge': 0, 'chk':'t1'}
    #ret = make_gaussian_input(data,fp_st)
    #with open('t1.com','w') as fp:
    #    fp.write(ret)

    lebesgue_context = LebesgueContext(
        executor="lebesgue_v2",
        extra='{"scass_type":"c16_m64_cpu","program_id":2118}',
        authorization='eyJhbGciOiJIUzI1NiIsInR5cCI6IkpXVCJ9.eyJleHAiOjE2NjIyODc3OTAsImlhdCI6MTY1OTY5NTc5MCwiaWRlbnRpdHkiOnsidXNlcl9pZCI6MTMyNCwidXNlcl9uYW1lIjoiMjM2OTU2MTMwMEBxcS5jb20iLCJlbWFpbCI6IjIzNjk1NjEzMDBAcXEuY29tIiwidXNlcl9raW5kIjoyLCJvcmdJZCI6MTMyNH19.mNFtXhaY7a4c-7E-vUQsaLy02nC6BmoqNfsXrL42R00',
        app_name='Default',
        org_id='123',
        user_id='456',
        tag='',
    )

    hello = ShellOPTemplate(name='Hello',
            image='LBG_Gaussian_1_v2',
            command=["bash"],
            script="cd /tmp/workdir; ulimit -s unlimited; source /root/g16.sh; g16 s0.com > s0.log")

    hello.inputs.artifacts = {"input": InputArtifact(path="/tmp/workdir")}
    helloArtiFact = upload_artifact(["./s0.com"])
    hello.outputs.artifacts = {"bar": OutputArtifact(path="/tmp/workdir/s0.log")}
    wf = Workflow(name="steps-lebesgue",context=lebesgue_context,host="http://39.106.93.187:32746")
    hello0 = Step(name="hello0",template=hello,artifacts={"input":helloArtiFact})
    wf.add(hello0)
    wf.submit()
    while wf.query_status() in ["Pending","Running"]:
        time.sleep(4)
    assert(wf.query_status() == 'Succeeded')
    step = wf.query_step(name="hello0")[0]
    download_artifact(step.outputs.artifacts["bar"])

    s0data = s0data('s0.log')
    np.save('s0data.npy',s0data)
    fp_s1 = {'nproc':16,'mem':4,'keywords':"opt b3lyp/gen pseudo=read TD",'pseudo':True,\
            'pseudo_element':'Ir','other_element':'C H O N','pseudo_basis':'Lanl2DZ',\
            'other_basis':'6-31G**','multiplicity': 1, 'charge': 0, 'chk':'s1'}
    ret = make_gaussian_input(s0data,fp_s1)
    with open('s1.com','w') as fp:
        fp.write(ret)
    fp_st = {'nproc':16,'mem':4,'keywords':"opt b3lyp/gen pseudo=read",'pseudo':True,\
            'pseudo_element':'Ir','other_element':'C H O N','pseudo_basis':'Lanl2DZ',\
            'other_basis':'6-31G**','multiplicity': 3, 'charge': 0, 'chk':'t1'}
    ret = make_gaussian_input(s0data,fp_st)
    with open('t1.com','w') as fp:
        fp.write(ret)
    hello = ShellOPTemplate(name='Hello',
            image='LBG_Gaussian_1_v2',
            command=["bash"],
            script="cd /tmp/workdir; ulimit -s unlimited; source /root/g16.sh; g16 s1.com > s1.log")

    hello.inputs.artifacts = {"input": InputArtifact(path="/tmp/workdir")}
    helloArtiFact = upload_artifact(["./s1.com"])
    hello.outputs.artifacts = {"bar": OutputArtifact(path="/tmp/workdir/s1.log")}
    wf = Workflow(name="steps-lebesgue",context=lebesgue_context,host="http://39.106.93.187:32746")
    hello0 = Step(name="hello0",template=hello,artifacts={"input":helloArtiFact})
    wf.add(hello0)
    wf.submit()
    while wf.query_status() in ["Pending","Running"]:
        time.sleep(4)
    assert(wf.query_status() == 'Succeeded')
    step = wf.query_step(name="hello0")[0]
    download_artifact(step.outputs.artifacts["bar"])

    hello = ShellOPTemplate(name='Hello',
            image='LBG_Gaussian_1_v2',
            command=["bash"],
            script="cd /tmp/workdir; ulimit -s unlimited; source /root/g16.sh; g16 t1.com > t1.log")

    hello.inputs.artifacts = {"input": InputArtifact(path="/tmp/workdir")}
    helloArtiFact = upload_artifact(["./t1.com"])
    hello.outputs.artifacts = {"bar": OutputArtifact(path="/tmp/workdir/t1.log")}
    wf = Workflow(name="steps-lebesgue",context=lebesgue_context,host="http://39.106.93.187:32746")
    hello0 = Step(name="hello0",template=hello,artifacts={"input":helloArtiFact})
    wf.add(hello0)
    wf.submit()
    while wf.query_status() in ["Pending","Running"]:
        time.sleep(4)
    assert(wf.query_status() == 'Succeeded')
    step = wf.query_step(name="hello0")[0]
    download_artifact(step.outputs.artifacts["bar"])

    s1data = s1data('s1.log')
    t1data = t1data('t1.log')
    np.save('s1data.npy',s1data)
    np.save('t1data.npy',t1data)
