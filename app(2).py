import streamlit as st
import pandas as pd
from io import BytesIO
from rdkit import Chem
from rdkit.Chem import Draw, Descriptors, Crippen, Lipinski, rdMolDescriptors, AllChem
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle
from reportlab.lib import colors
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.pagesizes import A4
import py3Dmol
import streamlit.components.v1 as components

st.set_page_config(page_title='ChemForge | Molecular Lab', page_icon='⚗️', layout='wide')
st.markdown('''<style>
.stApp {background:#0b1220;color:#e2e8f0}
[data-testid="stSidebar"] {background:#111c30}
h1,h2,h3 {color:#67e8f9!important}
div[data-testid="stMetric"] {background:#17253b;padding:15px;border-radius:14px;border:1px solid #28445b}
.stButton button {background:#0891b2;color:white;border-radius:10px;border:0}
</style>''',unsafe_allow_html=True)

EXAMPLES={'Ethanol':'CCO','Benzene':'c1ccccc1','Acetic acid':'CC(=O)O','Caffeine':'Cn1c(=O)c2c(ncn2C)n(C)c1=O','Aspirin':'CC(=O)Oc1ccccc1C(=O)O','Acetone':'CC(=O)C','Glycerol':'OCC(O)CO','Phenol':'Oc1ccccc1'}
if 'saved' not in st.session_state: st.session_state.saved={}
if 'smiles' not in st.session_state: st.session_state.smiles='CCO'

def props(m):
    h=Chem.AddHs(m)
    kek=Chem.Mol(m)
    try: Chem.Kekulize(kek,clearAromaticFlags=True)
    except Exception: pass
    pi=sum(2 if b.GetBondType()==Chem.BondType.TRIPLE else 1 if b.GetBondType()==Chem.BondType.DOUBLE else 0 for b in kek.GetBonds())
    return {'Formula':rdMolDescriptors.CalcMolFormula(m),'Molecular weight (g/mol)':round(Descriptors.MolWt(m),3),'LogP':round(Crippen.MolLogP(m),3),'TPSA (Å²)':round(rdMolDescriptors.CalcTPSA(m),2),'H-bond donors':Lipinski.NumHDonors(m),'H-bond acceptors':Lipinski.NumHAcceptors(m),'Rotatable bonds':Lipinski.NumRotatableBonds(m),'Rings':Lipinski.RingCount(m),'Heavy atoms':m.GetNumHeavyAtoms(),'Sigma bonds':h.GetNumBonds(),'Pi bonds':pi}

def image_bytes(m):
    out=BytesIO();Draw.MolToImage(m,size=(750,500)).save(out,format='PNG');return out.getvalue()

def report(name,smiles,p):
    out=BytesIO();doc=SimpleDocTemplate(out,pagesize=A4)
    styles=getSampleStyleSheet();rows=[['Descriptor','Value']]+[[k,str(v)] for k,v in p.items()]
    t=Table(rows,colWidths=[260,220],repeatRows=1)
    t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#0891b2')),('TEXTCOLOR',(0,0),(-1,0),colors.white),('GRID',(0,0),(-1,-1),0.4,colors.grey),('PADDING',(0,0),(-1,-1),8)]))
    doc.build([Paragraph('ChemForge | Molecular Analysis',styles['Title']),Spacer(1,20),Paragraph('Student: '+name,styles['Normal']),Paragraph('SMILES: '+smiles,styles['Normal']),Spacer(1,20),t,Spacer(1,20),Paragraph('Computed descriptors are educational predictions, not experimental measurements.',styles['Normal'])]);return out.getvalue()

def mol_from(s):
    try: return Chem.MolFromSmiles(s.strip()) if s.strip() else None
    except Exception: return None

st.sidebar.title('⚗️ ChemForge')
page=st.sidebar.radio('Workspace',['🏠 Dashboard','🧬 Molecule Explorer','📊 Compare Molecules','📚 My Collection','🧠 Knowledge Check'])
student=st.sidebar.text_input('Student name',value='Student')
st.sidebar.caption('Interactive molecular science • Powered by RDKit')

if page=='🏠 Dashboard':
    st.title('⚗️ ChemForge Molecular Laboratory')
    st.write('Explore structures, calculate molecular descriptors, and compare chemical compounds.')
    a,b,c=st.columns(3)
    a.metric('Example compounds',len(EXAMPLES));b.metric('Saved this session',len(st.session_state.saved));c.metric('Analysis modes','2D + 3D')
    st.subheader('Quick-start molecules')
    for name,s in EXAMPLES.items():
        if st.button(f'Explore {name} · {s}',key='home'+name):
            st.session_state.smiles=s;st.session_state.go_explorer=True;st.rerun()
    st.info('Choose Molecule Explorer from the sidebar to inspect a molecule. Your collection is stored only for this browser session.')

elif page=='🧬 Molecule Explorer':
    st.title('🧬 Molecule Explorer')
    ex=st.selectbox('Load an example', ['Custom']+list(EXAMPLES))
    if st.button('Load selected example') and ex!='Custom':st.session_state.smiles=EXAMPLES[ex]
    s=st.text_input('Enter SMILES notation',key='smiles')
    m=mol_from(s)
    if m is None: st.warning('Enter a valid SMILES string to continue.')
    else:
        p=props(m)
        st.success('Valid molecular structure')
        x,y,z=st.columns(3);x.metric('Formula',p['Formula']);y.metric('Mass (g/mol)',p['Molecular weight (g/mol)']);z.metric('LogP',p['LogP'])
        left,right=st.columns(2)
        with left:
            st.subheader('2D structure');st.image(image_bytes(m),use_container_width=True)
            st.download_button('Download 2D PNG',image_bytes(m),file_name='chemforge_structure.png',mime='image/png')
        with right:
            st.subheader('Interactive 3D structure')
            style=st.selectbox('Display style',['Ball and stick','Stick','Space filling','Wireframe'])
            m3=Chem.AddHs(m);params=AllChem.ETKDGv3();params.randomSeed=42
            if AllChem.EmbedMolecule(m3,params)==0:
                try: AllChem.MMFFOptimizeMolecule(m3,maxIters=150)
                except Exception: pass
                v=py3Dmol.view(width=650,height=450)
                v.addModel(Chem.MolToMolBlock(m3),'mol')
                options={'Ball and stick':{'stick':{'radius':0.15},'sphere':{'scale':0.28}},'Stick':{'stick':{}},'Space filling':{'sphere':{}},'Wireframe':{'line':{}}}
                v.setStyle(options[style]);v.setBackgroundColor('#ffffff');v.zoomTo()
                components.html(v._make_html(),height=470)
            else: st.info('3D coordinates could not be generated for this molecule.')
        st.subheader('Molecular descriptors')
        st.dataframe(pd.DataFrame([{'Descriptor':k,'Value':v} for k,v in p.items()]),hide_index=True,use_container_width=True)
        st.download_button('Download properties CSV',pd.DataFrame([p]).to_csv(index=False).encode(),file_name='chemforge_properties.csv',mime='text/csv')
        if st.button('⭐ Save molecule to collection'):
            st.session_state.saved[Chem.MolToSmiles(m)]=p;st.success('Saved to this session.')
        st.download_button('📄 Download PDF analysis',report(student,s,p),file_name='chemforge_report.pdf',mime='application/pdf')

elif page=='📊 Compare Molecules':
    st.title('📊 Compare Molecules')
    a,b=st.columns(2)
    s1=a.text_input('First SMILES',value='CCO');s2=b.text_input('Second SMILES',value='CC(=O)O')
    m1,m2=mol_from(s1),mol_from(s2)
    if m1 is None or m2 is None: st.error('Both SMILES must be valid.')
    else:
        p1,p2=props(m1),props(m2)
        a.image(image_bytes(m1),caption='Molecule A',use_container_width=True)
        b.image(image_bytes(m2),caption='Molecule B',use_container_width=True)
        df=pd.DataFrame({'Descriptor':list(p1),'Molecule A':list(p1.values()),'Molecule B':list(p2.values())})
        st.dataframe(df,hide_index=True,use_container_width=True)
        st.download_button('Download comparison CSV',df.to_csv(index=False).encode(),file_name='comparison.csv',mime='text/csv')
        st.caption('A descriptor difference does not alone establish safety, efficacy, or experimental behavior.')

elif page=='📚 My Collection':
    st.title('📚 My Molecule Collection')
    if not st.session_state.saved:st.info('Save molecules from the Explorer to see them here.')
    else:
        df=pd.DataFrame([{'SMILES':s,**p} for s,p in st.session_state.saved.items()])
        st.dataframe(df,hide_index=True,use_container_width=True)
        st.download_button('Export collection CSV',df.to_csv(index=False).encode(),file_name='my_molecules.csv',mime='text/csv')
        if st.button('Clear collection'):st.session_state.saved={};st.rerun()

else:
    st.title('🧠 Knowledge Check')
    questions=[('What does SMILES represent?',['A text notation for molecular structures','A microscope','A laboratory chemical'],'A text notation for molecular structures'),('Which descriptor estimates lipophilicity?',['TPSA','LogP','Ring count'],'LogP'),('Which library calculates molecular descriptors here?',['RDKit','Pygame','Tkinter'],'RDKit')]
    with st.form('quiz'):
        answers=[st.radio(q,options,key=f'q{i}') for i,(q,options,_) in enumerate(questions)]
        submitted=st.form_submit_button('Check answers')
    if submitted:
        score=sum(answers[i]==correct for i,(_,_,correct) in enumerate(questions))
        st.success(f'Score: {score}/{len(questions)}')
        for i,(q,_,correct) in enumerate(questions):
            if answers[i]!=correct:st.write(f'Question {i+1}: correct answer — {correct}')
