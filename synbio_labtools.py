import streamlit as st
# To make things easier later, we're also importing numpy and pandas for
# working with sample data.
import numpy as np
import pandas as pd
import base64


def get_table_download_link(df):
    """Generates a link allowing the data in a given panda dataframe to be downloaded
    in:  dataframe
    out: href string
    """
    csv = df.to_csv(index=True)
    b64 = base64.b64encode(csv.encode()).decode()  # some strings <-> bytes conversions necessary here
    href = f'<a href="data:file/csv;base64,{b64}" download="gibson.csv">Download csv file</a>'
    return href

st.title('Gibson Calculator')

insert_number = st.number_input("Number of inserts", key='start', value=1)
vector_mass = st.number_input("Desired vector mass to use", key='vector_mass', value = 30.00)

length = dict()
concentration = dict()
ratio = dict()
name = dict()

column1, column2, column3 = st.beta_columns(3)

name['vector'] = column1.text_input('Vector', key='vector', value='vector')
length['vector'] = column2.number_input('Lenght (bp)', key='999', value=1000)
concentration['vector'] = column3.number_input('Concentration (ng/ul)', key='999', value=30.0)

ratio['vector'] = 1

c1, c2, c3, c4 = st.beta_columns(4)

for i in range(1, insert_number):
    name['insert'+str(i)] = c1.text_input(f'Insert {str(i)}', key=str(i), value='insert '+str(i))
    length['insert'+str(i)] = c2.number_input('Lenght (bp)', key=str(i), value=100)
    concentration['insert'+str(i)] = c3.number_input('Concentration (ng/ul)', key='lol'+str(i), value=1.0)
    ratio['insert'+str(i)] = c4.number_input('Insert:vector ratio', key='lol'+str(i), value=3)



data = pd.DataFrame([name, length, concentration, ratio]).T
data.columns = ['name','length' , 'conc', 'ratio']

data['mass'] = vector_mass * data['ratio'] * (data['length']/max(data['length']))
data['vol'] = data['mass']/data['conc']

data.set_index('name', drop=True, inplace=True)

data.columns = ['Length (bp)', "Concentration (ng/ul)",
                  "insert:vector ratio", "mass (ng)", "volume (ul)"]

data = data.iloc[:,[3,4,0,1,2]]

st.table(data)

st.markdown(get_table_download_link(data), unsafe_allow_html=True)



