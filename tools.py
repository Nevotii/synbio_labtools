# -*- coding: utf-8 -*-

import streamlit as st
import numpy as np
import pandas as pd
import base64
import fuzzysearch as fz
from Bio import SeqIO
from Bio.Seq import Seq


def get_table_download_link(df):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=True)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}" download="gibson.csv">Download table (.csv)</a>'
    return href

def intro():
    st.sidebar.success("Select a tool.")

    st.markdown(
        """
        Here is a set of tools created to ease the life of a synbio lab
        
        **👈 Select a tool from the dropdown on the left**
    """
    )

def gibson():

    st.title('Gibson Calculator')

    insert_number = st.number_input("Number of parts", key='start', value=1, min_value=1, step=1)
    vector_mass = st.number_input("Desired vector mass to use (ng)", key='vector_mass',min_value=0.0, value = 30.00)

    length = dict()
    concentration = dict()
    ratio = dict()
    name = dict()

    column1, column2, column3 = st.columns(3)

    name['vector'] = column1.text_input('Vector', key='vector', value='vector')
    length['vector'] = column2.number_input('Lenght (bp)', key=''vectorvector, value=1000, min_value=0)
    concentration['vector'] = column3.number_input('Concentration (ng/ul)', key='concconc', value=30.0, min_value=0.0)

    ratio['vector'] = 1

    c1, c2, c3, c4 = st.columns(4)

    for i in range(1, insert_number):
        name['insert'+str(i)] = c1.text_input(f'Insert {str(i)}', key=str(i), value='insert '+str(i))
        length['insert'+str(i)] = c2.number_input('Lenght (bp)', key=str(i), value=100, min_value=0)
        concentration['insert'+str(i)] = c3.number_input('Concentration (ng/ul)', key='lol'+str(i), value=1.0, min_value=0.0)
        ratio['insert'+str(i)] = c4.number_input('Insert:vector ratio', key='lol'+str(i), value=3.0, min_value=0.0)



    data = pd.DataFrame([name, length, concentration, ratio]).T
    data.columns = ['name','length' , 'conc', 'ratio']

    data['mass'] = vector_mass * data['ratio'] * (data['length']/max(data['length']))
    data['vol'] = data['mass']/data['conc']
    data['vol'] = data['vol'].astype('float').round(2)
    total_volume = data['vol'].sum()
    

    data.set_index('name', drop=True, inplace=True)
    data.loc['Total'] = pd.Series(data['vol'].sum(), index = ['vol'])
    data = data.fillna('-')
    
    data.columns = ['Length (bp)', "Concentration (ng/ul)",
                      "insert:vector ratio", "mass (ng)", "volume (ul)"]

    data = data.iloc[:,[3,4,0,1,2]]
    

    st.table(data)
    st.markdown(get_table_download_link(data), unsafe_allow_html=True)


    st.subheader('References')
    st.write('Mass of inserts is calculated using the following formula:')
    st.latex(r'mass\ of\ insert\ =\ insert:vector \ molar \ ratio * mass\ of\ vector * \frac{insert \ length}{vector \ length}')



def pattern():
    with open('cutibacterium_acnes.fasta', 'r') as genome:
        seq = SeqIO.read(genome, 'fasta').seq

    sequence = st.text_area('Sequence to search in C. acnes genome')
    mismatch = st.number_input('Select number of mismatches', value=0, min_value=0, max_value=5)

    if len(sequence) < 15:
        st.error('Please write a sequence longer than 15bp')
    subseq = Seq(sequence.upper())

    if sequence:
        fw = fz.find_near_matches(subseq, seq, max_l_dist=mismatch, max_deletions=0,max_insertions=0)
        rv = fz.find_near_matches(subseq.reverse_complement(), seq, max_l_dist=mismatch, max_deletions=0,max_insertions=0)
        
        if fw or rv:
            st.write(f'## You have {len(fw)+len(rv)} matches for your sequence in KPA')
            c1, c2, c3, c4 = st.columns(4)
            c1.write('### Position')
            c2.write('### Mismatch')
            c3.write('### Sequence')
            c4.write('### Direction')
            for match in fw:
                c1.write(match.start)
                c2.write(match.dist)
                c3.write(match.matched)
                c4.write('`forward`')
            for match in rv:
                c1.write(match.start)
                c2.write(match.dist)
                c3.write(match.matched)
                c4.write('`reverse`')
        else:
            st.error('No matches found')
            
def baseedit():
    """
    Created on Tue Jul 16 13:51:30 2019
    This programm shall find all spacer designs to insert intendet mutation into your genome
    @author: chdavo
    """
    
    import base_edit
    import io
    
    uploaded_file = st.file_uploader("Choose a file", key="genome")
    
    if uploaded_file is not None:
        genome = io.StringIO(uploaded_file.getvalue().decode("utf-8"))
        seq = SeqIO.read(genome, 'fasta').seq
    else:
        with open('cutibacterium_acnes.fasta', 'r') as genome:
            seq = SeqIO.read(genome, 'fasta').seq

    uploaded_fasta = st.file_uploader("Choose a file ORF", key="ORF")    
    
    if uploaded_fasta is not None:
        fasta = io.StringIO(uploaded_fasta.getvalue().decode("utf-8"))
        count_targetedgenes = 0
        placestopcodon=[]
        print("I start iterating...")
        for seq_record in SeqIO.parse(fasta, 'fasta'):
            print(seq_record.seq)
            placestopcodon = base_edit.Spacer_search(seq_record,placestopcodon)[1]
            print(placestopcodon)
            
        for item in placestopcodon:
            st.write("%s\n" % item)
    else:
        st.write("funcion acabada")
    