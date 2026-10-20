#!/usr/bin/env python
# coding: utf-8

# In[50]:


import pandas as pd
import scipy.stats as stats
import tabulate as tabulate
from tabulate import tabulate
import numpy as np
from scipy.stats import f_oneway, ttest_ind, chi2_contingency
from statsmodels.stats.multicomp import pairwise_tukeyhsd
import regex as re

#file_xlsx = './data/ADNIMERGE_230712.xlsx'
file_xlsx = './data/ADNIMERGE_220706.xlsx'
file_json = './data/categories_nm.json'


# In[51]:


df=pd.read_excel(file_xlsx)


# In[52]:


#NM subjects
data_nm={}
data_nm["Group"]=["Subjects", "Visits", "Female", "APOE4", "Age", "Education", "MMSE", "CDRSB", "ADAS13", "FAQ", "LDELTOTAL", "TRABSCORE", "PACC", "RAVLT-Im", "RAVLT-Lear", "RAVLT-perFor"]

for dx_category in ["CN", "SMC", "MCI", "AD",]:
    #First filtering (by group)
    if dx_category == "MCI":
        filtered_df = df[((df["DX_bl"] == "EMCI") | (df["DX_bl"] == "LMCI"))]
    else:
        filtered_df = df[(df["DX_bl"] == dx_category)]
   
    #VISITS
    visits = filtered_df.dropna(subset="DX").shape[0]
    visits_percentage = f"{(visits/df.dropna(subset='DX').shape[0])*100:.1f}%"
   
    #Second filtering (at baseline)
    filtered_df = filtered_df[filtered_df["VISCODE"] == "bl"]
   
    #SUBJECTS
    unique_subjects = filtered_df["RID"].nunique()
    subjects = df["RID"].nunique()
    unique_subjects_percentage = f"{(unique_subjects/subjects)*100:.1f}%"
    
    #FEMALE
    female_subjects = filtered_df[filtered_df["PTGENDER"] == "Female"]["RID"].nunique()
    female_percentage = f"{(female_subjects/unique_subjects)*100:.1f}%"
   
    #APOE4
    apoE4_subjects = filtered_df[filtered_df["APOE4"] > 0]["RID"].nunique()
    apoE4_n = filtered_df[(filtered_df["APOE4"] >=0) & (filtered_df["APOE4"]<=2)]["RID"].nunique()
    apoE4_percentage = f"{(apoE4_subjects/apoE4_n)*100:.1f}%"
    
    #AGE
    age_mean = filtered_df["AGE"].mean()
    age_std = filtered_df["AGE"].std()
    age_min = filtered_df["AGE"].min()
    age_max = filtered_df["AGE"].max()
    age_stats = f"{age_mean:.1f} ({age_std:.1f}) ({age_min:.1f}, {age_max:.1f})"
   
    #EDUCATION
    education_mean = filtered_df["PTEDUCAT"].mean()
    education_std = filtered_df["PTEDUCAT"].std()
    education_min = filtered_df["PTEDUCAT"].min()
    education_max = filtered_df["PTEDUCAT"].max()
    education_stats = f"{education_mean:.1f} ({education_std:.1f}) ({education_min:.1f}, {education_max:.1f})"
   
    #MMSE
    mmse_mean = filtered_df["MMSE_bl"].mean()
    mmse_std = filtered_df["MMSE_bl"].std()
    mmse_min = filtered_df["MMSE_bl"].min()
    mmse_max = filtered_df["MMSE_bl"].max()
    mmse_stats = f"{mmse_mean:.1f} ({mmse_std:.1f}) ({mmse_min:.1f}, {mmse_max:.1f})"
  
    #CDRSB
    cdrsb_mean = filtered_df["CDRSB_bl"].mean()
    cdrsb_std = filtered_df["CDRSB_bl"].std()
    cdrsb_min = filtered_df["CDRSB_bl"].min()
    cdrsb_max = filtered_df["CDRSB_bl"].max()
    cdrsb_stats = f"{cdrsb_mean:.2f} ({cdrsb_std:.2f}) ({cdrsb_min:.2f}, {cdrsb_max:.2f})"
   
    #ADAS13
    adas13_mean = filtered_df["ADAS13_bl"].mean()
    adas13_std = filtered_df["ADAS13_bl"].std()
    adas13_min = filtered_df["ADAS13_bl"].min()
    adas13_max = filtered_df["ADAS13_bl"].max()
    adas13_stats = f"{adas13_mean:.1f} ({adas13_std:.1f}) ({adas13_min:.1f}, {adas13_max:.1f})"
   
    #FAQ
    faq_mean = filtered_df["FAQ_bl"].mean()
    faq_std = filtered_df["FAQ_bl"].std()
    faq_min = filtered_df["FAQ_bl"].min()
    faq_max = filtered_df["FAQ_bl"].max()
    faq_stats = f"{faq_mean:.1f} ({faq_std:.1f}) ({faq_min:.1f}, {faq_max:.1f})"
   
    #LDELTOTAL
    ldeltotal_mean = filtered_df["LDELTOTAL_BL"].mean()
    ldeltotal_std = filtered_df["LDELTOTAL_BL"].std()
    ldeltotal_min = filtered_df["LDELTOTAL_BL"].min()
    ldeltotal_max = filtered_df["LDELTOTAL_BL"].max()
    ldeltotal_stats = f"{ldeltotal_mean:.1f} ({ldeltotal_std:.1f}) ({ldeltotal_min:.1f}, {ldeltotal_max:.1f})"
    
    #TRABSCORE
    trabscore_mean = filtered_df["TRABSCOR_bl"].mean()
    trabscore_std = filtered_df["TRABSCOR_bl"].std()
    trabscore_min = filtered_df["TRABSCOR_bl"].min()
    trabscore_max = filtered_df["TRABSCOR_bl"].max()
    trabscore_stats = f"{trabscore_mean:.1f} ({trabscore_std:.1f}) ({trabscore_min:.1f}, {trabscore_max:.1f})"
    
    #PACC (mPACCdigit)
    pacc_mean = filtered_df["mPACCdigit_bl"].mean()
    pacc_std = filtered_df["mPACCdigit_bl"].std()
    pacc_min = filtered_df["mPACCdigit_bl"].min()
    pacc_max = filtered_df["mPACCdigit_bl"].max()
    pacc_stats = f"{pacc_mean:.1f} ({pacc_std:.1f}) ({pacc_min:.1f}, {pacc_max:.1f})"
   
    #RAVLT-immediate
    ravlt_im_mean = filtered_df["RAVLT_immediate_bl"].mean()
    ravlt_im_std = filtered_df["RAVLT_immediate_bl"].std()
    ravlt_im_min = filtered_df["RAVLT_immediate_bl"].min()
    ravlt_im_max = filtered_df["RAVLT_immediate_bl"].max()
    ravlt_im_stats = f"{ravlt_im_mean:.1f} ({ravlt_im_std:.1f}) ({ravlt_im_min:.1f}, {ravlt_im_max:.1f})"
   
    #RAVLT-learning
    ravlt_lear_mean = filtered_df["RAVLT_learning_bl"].mean()
    ravlt_lear_std = filtered_df["RAVLT_learning_bl"].std()
    ravlt_lear_min = filtered_df["RAVLT_learning_bl"].min()
    ravlt_lear_max = filtered_df["RAVLT_learning_bl"].max()
    ravlt_lear_stats = f"{ravlt_lear_mean:.1f} ({ravlt_lear_std:.1f}) ({ravlt_lear_min:.1f}, {ravlt_lear_max:.1f})"
   
    #RAVLT-percentage forgetting
    ravlt_perfor_mean = filtered_df["RAVLT_perc_forgetting_bl"].mean()
    ravlt_perfor_std = filtered_df["RAVLT_perc_forgetting_bl"].std()
    ravlt_perfor_min = filtered_df["RAVLT_perc_forgetting_bl"].min()
    ravlt_perfor_max = filtered_df["RAVLT_perc_forgetting_bl"].max()
    ravlt_perfor_stats = f"{ravlt_perfor_mean:.1f} ({ravlt_perfor_std:.1f}) ({ravlt_perfor_min:.1f}, {ravlt_perfor_max:.1f})"

    data_nm[dx_category]=[f"{unique_subjects} ({unique_subjects_percentage})", f"{visits} ({visits_percentage})", f"{female_subjects} ({female_percentage})", f"{apoE4_subjects} (n={apoE4_n} {apoE4_percentage})", age_stats, education_stats, mmse_stats, cdrsb_stats, adas13_stats, faq_stats, ldeltotal_stats, trabscore_stats, pacc_stats, ravlt_im_stats, ravlt_lear_stats, ravlt_perfor_stats]
    


# In[53]:


#Group1 : CN, Group2 : SMC, Group3 : MCI, Group4 : AD
#CN to SMC (a)
#CN to MCI (b)
#CN to AD (c)
#SMC to MCI (d)
#SMC to AD (e)
#MCI to AD (f)
def one_way_anova_with_posthoc(group1, group2, group3, group4):
    # Perform one-way ANOVA
    f_value, p_value = f_oneway(group1, group2, group3, group4)
    print(p_value)
    # Check if p-value is less than 0.01
    if p_value < 0.01:
        # Perform Tukey post hoc test
        data = np.concatenate([group1, group2, group3, group4])
        labels = ['group1'] * len(group1) + ['group2'] * len(group2) + ['group3'] * len(group3) + ['group4'] * len(group4)
        tukey_results = pairwise_tukeyhsd(data, labels, 0.01)
        possible_results = ['a','b','c','d','e','f']
        post_hoc = [res for res,reject in zip(possible_results,tukey_results.reject) if reject == True]
        post_hoc = ",".join(post_hoc)
        return post_hoc


# In[54]:


#Group1 : CN, Group2 : SMC, Group3 : MCI, Group4 : AD
#CN to SMC (a)
#CN to MCI (b)
#CN to AD (c)
#SMC to MCI (d)
#SMC to AD (e)
#MCI to AD (f)
def chi_square_with_posthoc(group1, group2, group3, group4, contingency_table):
    # Perform chi-square test
    chi2, p, dof, expected = chi2_contingency(contingency_table)

    # Perform post-hoc test if p-value is significant
    if p < 0.01:
        data = np.concatenate([group1, group2, group3, group4])
        labels = ['group1'] * len(group1) + ['group2'] * len(group2) + ['group3'] * len(group3) + ['group4'] * len(group4)
        tukey_results = pairwise_tukeyhsd(data, labels, 0.01)
        possible_results = ['a','b','c','d','e','f']
        post_hoc = [res for res,reject in zip(possible_results,tukey_results.reject) if reject == True]
        post_hoc = ",".join(post_hoc)
        return post_hoc


# In[55]:


filtered_nm_CN = df[(df["DX_bl"] == "CN") & (df["VISCODE"] == "bl")]
filtered_nm_SMC = df[(df["DX_bl"] == "SMC") & (df["VISCODE"] == "bl")]
filtered_nm_MCI = df[((df["DX_bl"] == "LMCI") | (df["DX_bl"] == "EMCI")) & (df["VISCODE"] == "bl")]
filtered_nm_AD = df[(df["DX_bl"] == "AD") & (df["VISCODE"] == "bl")]

data_nm["Post-Hoc"] = [""]*len(data_nm["Group"])


# In[56]:


group1_nm = filtered_nm_CN["PTGENDER"].dropna().replace({'Male': 1, 'Female': 0})
group2_nm = filtered_nm_SMC["PTGENDER"].dropna().replace({'Male': 1, 'Female': 0})
group3_nm = filtered_nm_MCI["PTGENDER"].dropna().replace({'Male': 1, 'Female': 0})
group4_nm = filtered_nm_AD["PTGENDER"].dropna().replace({'Male': 1, 'Female': 0})


# Convert data to contingency table format
contingency_table_nm = np.array([
    [np.sum(group1_nm == 0), np.sum(group1_nm == 1)],
    [np.sum(group2_nm == 0), np.sum(group2_nm == 1)],
    [np.sum(group3_nm == 0), np.sum(group3_nm == 1)],
    [np.sum(group4_nm == 0), np.sum(group4_nm == 1)]
])

data_nm["Post-Hoc"][data_nm["Group"].index("Female")] = chi_square_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm, contingency_table_nm)


# In[57]:


#APOE4
group1_nm = filtered_nm_CN["APOE4"].dropna()
group2_nm = filtered_nm_SMC["APOE4"].dropna()
group3_nm = filtered_nm_MCI["APOE4"].dropna()
group4_nm = filtered_nm_AD["APOE4"].dropna()

# Convert data to contingency table format
contingency_table_nm = np.array([
    [np.sum(group1_nm == 0), np.sum(group1_nm == 1), np.sum(group1_nm == 2)],
    [np.sum(group2_nm == 0), np.sum(group2_nm == 1), np.sum(group1_nm == 2)],
    [np.sum(group3_nm == 0), np.sum(group3_nm == 1), np.sum(group1_nm == 2)],
    [np.sum(group4_nm == 0), np.sum(group4_nm == 1), np.sum(group1_nm == 2)]
])

data_nm["Post-Hoc"][data_nm["Group"].index("APOE4")] = chi_square_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm, contingency_table_nm)


# In[58]:


#AGE
group1_nm = filtered_nm_CN["AGE"].dropna()
group2_nm = filtered_nm_SMC["AGE"].dropna()
group3_nm = filtered_nm_MCI["AGE"].dropna()
group4_nm = filtered_nm_AD["AGE"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("Age")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[59]:


#EDUCATION
group1_nm = filtered_nm_CN["PTEDUCAT"].dropna()
group2_nm = filtered_nm_SMC["PTEDUCAT"].dropna()
group3_nm = filtered_nm_MCI["PTEDUCAT"].dropna()
group4_nm = filtered_nm_AD["PTEDUCAT"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("Education")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[60]:


#MMSE
group1_nm = filtered_nm_CN["MMSE_bl"].dropna()
group2_nm = filtered_nm_SMC["MMSE_bl"].dropna()
group3_nm = filtered_nm_MCI["MMSE_bl"].dropna()
group4_nm = filtered_nm_AD["MMSE_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("MMSE")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[61]:


#CDRSB
group1_nm = filtered_nm_CN["CDRSB_bl"].dropna()
group2_nm = filtered_nm_SMC["CDRSB_bl"].dropna()
group3_nm = filtered_nm_MCI["CDRSB_bl"].dropna()
group4_nm = filtered_nm_AD["CDRSB_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("CDRSB")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[62]:


#ADAS13
group1_nm = filtered_nm_CN["ADAS13_bl"].dropna()
group2_nm = filtered_nm_SMC["ADAS13_bl"].dropna()
group3_nm = filtered_nm_MCI["ADAS13_bl"].dropna()
group4_nm = filtered_nm_AD["ADAS13_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("ADAS13")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[63]:


#FAQ
group1_nm = filtered_nm_CN["FAQ_bl"].dropna()
group2_nm = filtered_nm_SMC["FAQ_bl"].dropna()
group3_nm = filtered_nm_MCI["FAQ_bl"].dropna()
group4_nm = filtered_nm_AD["FAQ_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("FAQ")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[64]:


#LDELTOTAL
group1_nm = filtered_nm_CN["LDELTOTAL_BL"].dropna()
group2_nm = filtered_nm_SMC["LDELTOTAL_BL"].dropna()
group3_nm = filtered_nm_MCI["LDELTOTAL_BL"].dropna()
group4_nm = filtered_nm_AD["LDELTOTAL_BL"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("LDELTOTAL")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[65]:


#TRABSCORE
group1_nm = filtered_nm_CN["TRABSCOR_bl"].dropna()
group2_nm = filtered_nm_SMC["TRABSCOR_bl"].dropna()
group3_nm = filtered_nm_MCI["TRABSCOR_bl"].dropna()
group4_nm = filtered_nm_AD["TRABSCOR_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("TRABSCORE")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[66]:


#PACC
group1_nm = filtered_nm_CN["mPACCdigit_bl"].dropna()
group2_nm = filtered_nm_SMC["mPACCdigit_bl"].dropna()
group3_nm = filtered_nm_MCI["mPACCdigit_bl"].dropna()
group4_nm = filtered_nm_AD["mPACCdigit_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("PACC")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[67]:


#RAVLT-Im
group1_nm = filtered_nm_CN["RAVLT_immediate_bl"].dropna()
group2_nm = filtered_nm_SMC["RAVLT_immediate_bl"].dropna()
group3_nm = filtered_nm_MCI["RAVLT_immediate_bl"].dropna()
group4_nm = filtered_nm_AD["RAVLT_immediate_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("RAVLT-Im")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[68]:


#RAVLT-Lear
group1_nm = filtered_nm_CN["RAVLT_learning_bl"].dropna()
group2_nm = filtered_nm_SMC["RAVLT_learning_bl"].dropna()
group3_nm = filtered_nm_MCI["RAVLT_learning_bl"].dropna()
group4_nm = filtered_nm_AD["RAVLT_learning_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("RAVLT-Lear")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[69]:


#RAVLT-perfor
group1_nm = filtered_nm_CN["RAVLT_perc_forgetting_bl"].dropna()
group2_nm = filtered_nm_SMC["RAVLT_perc_forgetting_bl"].dropna()
group3_nm = filtered_nm_MCI["RAVLT_perc_forgetting_bl"].dropna()
group4_nm = filtered_nm_AD["RAVLT_perc_forgetting_bl"].dropna()

data_nm["Post-Hoc"][data_nm["Group"].index("RAVLT-perFor")] = one_way_anova_with_posthoc(group1_nm, group2_nm, group3_nm, group4_nm)


# In[70]:


#Display table (and latex format)
data_nm["Post-Hoc"] = ["" if res == None else res for res in data_nm["Post-Hoc"]]
table_nm = pd.DataFrame(data_nm)
latex_table = table_nm.to_latex(index=False)
print(latex_table)
table_nm = table_nm.style.hide(axis="index")
display(table_nm)

