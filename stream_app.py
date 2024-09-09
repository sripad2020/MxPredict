import streamlit as st
import pandas as pd
import joblib
from collections import Counter
st.title('Monkey Pox Clade Prediction Tool')
def neuclotide_features(seq):
    cleaned_seq = ''.join([base for base in seq if base in 'ATCG'])
    length = len(cleaned_seq)
    a_count = cleaned_seq.count('A')
    t_count = cleaned_seq.count('T')
    c_count = cleaned_seq.count('C')
    g_count = cleaned_seq.count('G')
    gc_content = (g_count + c_count) / length * 100 if length > 0 else 0
    dinucleotide_count = Counter(cleaned_seq[i:i + 2] for i in range(len(cleaned_seq) - 1))
    trimer_count = Counter(cleaned_seq[i:i + 3] for i in range(len(cleaned_seq) - 2))
    features = {
        'length': length,
        'A_count': a_count,
        'T_count': t_count,
        'C_count': c_count,
        'G_count': g_count,
        'GC_content': gc_content,
    }
    for dinucleotide in ['AA', 'AT', 'TA', 'TT', 'TC', 'CT', 'CA', 'TG', 'GG', 'GA', 'GT', 'AC', 'AG', 'CC', 'CG',
                         'GC']:
        features[dinucleotide] = dinucleotide_count[dinucleotide]

    for trimer in ['AAA', 'AAT', 'ATA', 'TAT', 'ATT', 'TTT', 'TTC', 'TCT', 'CTT', 'CTA', 'ATC',
                   'TCA', 'CAT', 'ATG', 'TGG', 'GGA', 'GAT', 'TAA', 'TTG', 'TGT', 'GTG', 'GTC',
                   'CTG', 'GTA', 'TAC', 'ACA', 'AAG', 'AGT', 'GTT', 'TTA', 'TAG', 'ACC', 'CCA',
                   'CAC', 'ACT', 'TCG', 'CGC', 'GCT', 'TGC', 'GCG', 'CGT', 'GCC', 'CCC', 'CCT',
                   'CTC', 'CAA', 'AAC', 'CGA', 'CAG', 'AGC', 'GCA', 'ACG', 'AGA', 'AGG', 'TCC',
                   'GAG', 'GAA', 'TGA', 'GAC', 'CGG', 'GGC', 'GGT', 'CCG', 'GGG']:
        features[trimer] = trimer_count[trimer]
    return features

atgc=st.text_input("Enter the ATGC sequence")

if st.button("click for prediction") and atgc:
    col = atgc.split(',')
    col = pd.DataFrame(data=col, columns=['sequence'])
    neu_features = col['sequence'].apply(neuclotide_features).apply(pd.Series)
    print(len(neu_features.columns))
    neu_features = neu_features.drop(['GC_content'], axis=1)
    cols = neu_features.values[0]
    print(cols)
    model = joblib.load('Ext_predict.pkl')
    predictions = model.predict([cols])
    print(predictions)
    if predictions[0] == 0:
        result = "Its Clade Ia"
        st.write(result)
    elif predictions[0] == 1:
        result = "Its Clade Ib"
        st.write(result)
    elif predictions[0] == 2:
        result = "Its clade IIa"
        st.write(result)
    elif predictions[0] == 3:
        result = "Its clade IIb"
        st.write(result)