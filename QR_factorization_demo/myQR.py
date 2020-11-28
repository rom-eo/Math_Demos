#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:35:37 2020

@author: a
"""
#import sympy as sp#
from sympy import init_printing,latex,Matrix
init_printing()    
#
from IPython.display import display, Math

def texPrint(strg):
    display(Math(r"\text{"+strg+"}"))
    
def myQR(mat):
    
    
    m=mat.cols
    # displays the v
    temp=[]
    texPrint(r'''Problem. Find an orthonormal basis of the span of the vectors:''')
    for i in range(m): 
        temp.append(fr"v_{i+1} =")
        temp.append(latex(mat[:,i]))
        temp.append("; ")
    display(Math(" ".join(temp)))
        
    if len(mat.rref()[1])!=m:
        display(mat)
        texPrint("Matrix has dependent columns.")
    else:    
        w1=mat[:,0]
        w_norm=[]
        w_norm.append(w1.norm())
        u1=w1/w1.norm()
        u=[u1]
        texPrint("Solution. Set ")#,end="")
        temp=[fr" w_1=",latex(w1) ,';']
        display(Math(" ".join(temp)))
        texPrint("Then ")#,end='')
        temp=["\|w_1\| = ",latex(w1.norm()),r"; \ u_1:=\frac{1}{\|w_1\|}w_{1}=",latex(u1),'.']
        display(Math(" ".join(temp)))
    
        for i in range(1,m):
            temp=[]
            temp2=[]
            v_i=mat[:,i]
            temp.append(r"\text{Next, }")
            temp.append(fr"w_{i+1}=v_{i+1}")
            temp2.extend([fr"w_{i+1}= ",latex(v_i)])
            w_i=v_i
            for k,vec in enumerate(u):
                coef=v_i.dot(vec)
                temp.append("-")
                temp2.append("-")
                temp.append(fr"\langle v_{i+1},u_{k+1}\rangle u_{k+1}")
                temp2.extend([r"\left(",latex(coef),r"\right)"])
                temp2.append(latex(vec))
                w_i-=coef*vec
            #temp2.extend([" = ",latex(w_i)])
            display(Math(" ".join(temp)))
            display(Math(" ".join(temp2)))
            u_i=w_i/w_i.norm()
            w_norm.append(w_i.norm())
            temp=[fr"w_{i+1}=",latex(w_i),fr";\ \|w_{i+1}\| = ",latex(w_i.norm()),
                  fr";\ u_{i+1}=",latex(u_i)]
            display(Math(" ".join(temp)))
    
            u.append(u_i)
        display(Math(r"\text{An orthonormal basis is given by:}"))
        display(Math(latex(u)))
        # compute R
        R=Matrix.zeros(m)
        for j in range(m):
            R[j,j]=w_norm[j]
            for i in range(0,j):
                R[i,j]=mat[:,j].dot(u[i])
        # compute R repr   
        Rf=[['0' for _ in range(m)] for _ in range(m)]
        for j in range(m):
            Rf[j][j]=fr"\|w_{j+1}\|"
            for i in range(0,j):
                Rf[i][j]=fr"v_{j+1} \cdot u_{i+1}"#mat[:,j].dot(u[i])
        RfTex=r"\begin{bmatrix}"
        RfTex+=r"\\".join(["&".join(a) for a in Rf])
        RfTex+=r"\end{bmatrix}."
        
        texPrint("The Matrix R is")
        display(Math(RfTex))
        texPrint("That is")
        display(Math(latex(R)))
        U=u[0]
        for i in range(1,m):
            U=U.row_join(u[i])
        texPrint("The QR factorization  is")
        temp=[latex(mat),'=',latex(U),r'\cdot',latex(R)]
        display(Math(" ".join(temp)))
        #display(Math(latex(U*R-mat)))

def QRrandomMatrix():
    '''Pick a random Matrix
    and call myQR on it.'''
    import random
    n=random.choice(range(2,4))
    m=random.choice(range(2,4))
    pool=range(-5,5)
    mat=Matrix(n,m,[random.choice(pool) for _ in range(n*m)])
    myQR(mat)
    
if __name__=='__main__':    
    myQR(Matrix([[1,1,1,1],[0,1,2,1],[-1,-1,1,1]]).transpose())
    #myQR(Matrix([[1,1,0,0],[1,0,2,0],[2,2,0,3]]).transpose())