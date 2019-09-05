# -*- coding: utf-8 -*-
"""
Created on Fri May 24 10:32:58 2019

@author: fuyin

plot comparison 

"""

import numpy as np
import matplotlib.pyplot as plt


def read_file(fname):
    with open(fname, 'r') as f:
        x = f.readlines()
        l = [[float(word) for word in ll.split('\t')] for ll in x]
        return np.array(l)
    # return np.array(la)

def plot_comp(l1, l2):
    x = l1[:,1]
    et1 = l1[:, 2]
    et2 = l2[:, 2]
    plt.close('all')
    fig, ax1 = plt.subplots(1,1)
    fig.subplots_adjust(hspace=.7)
    ax1.plot(x, et1,'b-o', label='Sampson error')
    ax1.set_ylabel('translation error [%]')
    ax1.plot(x, et2, 'r-s', label='transfer error')
    
    # ax1.errorbar(x, m_ed, yerr = m_sd, fmt = 'ro')
    
    # ax1.plot(x, mean_da[:,1],'g-', label='')
    # ax1.plot(x, bias_a[:,2],'r-', label='bias_az')
    # ax1.set_ylim([0,40])
    # plt.title('Accelerometer Bias')
    ax1.set_xlabel('Feature point noise')
    ax1.set_title('Translation comparison transfer error vs Sampson error')
    ax1.legend(loc='upper right', fontsize='small')
    ax1.grid()
    plt.savefig("../result/translation_comparison.png", dpi=360)
    
    plt.draw()

if __name__=="__main__":
    l1 = read_file('../result/sampson_c_err.log')
    l2 = read_file('../result/transfer_err.log')
    plot_comp(l1, l2)
    # generate_data_for_test()
    # test_one_pose()
    
    
    
    
    
    
    
    