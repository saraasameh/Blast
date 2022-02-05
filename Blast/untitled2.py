# -*- coding: utf-8 -*-
"""
Created on Tue Dec 29 21:19:14 2020

@author: pc
"""

def matrix(a,b):
    with open("blosum62.txt") as matrix_file:
        matrix = matrix_file.read()
        lines = matrix.strip().split('\n')
    
        header = lines.pop(0)
        columns = header.split()
        matrix = {}
    
        for row in lines:
            entries = row.split()
            row_name = entries.pop(0)
            matrix[row_name] = {}
    
            if len(entries) != len(columns):
                raise Exception('Improper entry number in row')
            for column_name in columns:
                matrix[row_name][column_name] = entries.pop(0)
        
    a = a.upper()
    b = b.upper()

    if a not in matrix or b not in matrix[a]:
        raise Exception('[%s, %s]' % (a, b))
    return matrix[a][b]
#************************************************************
def remove_repeats(seq):
    NOTremoved=""
    index=0
    cnt=0
    while True:
        if (index+3 <= len(seq)):
            x=matrix( seq[index], seq[index])
            y=matrix( seq[index+1], seq[index+1])
            if matrix( seq[index+2], seq[index+2])==x and matrix( seq[index+3], seq[index+3])==y:
                while True:
                    if matrix( seq[index+2], seq[index+2])==matrix( seq[index], seq[index]) or matrix( seq[index+3], seq[index+3])==matrix( seq[index+1], seq[index+1]):
                        index+=2
                        cnt+=2
                    else:
                        index+=2
                        cnt+=2
                        if cnt==4:
                            for i in range (4): #remove if more than 4
                                NOTremoved+=seq[index-4+i]
                        cnt=0
                        break
                    if(index+2 == len(seq)): #if the next two are the last
                        index+=2
                        cnt+=2
                        if cnt==4:
                            for i in range (4):
                                NOTremoved+=seq[index-4+i]
                        cnt=0
                        break
            else:
                NOTremoved+=seq[index]
                index+=1
                
        elif (index == len(seq)):
            break
        else:
            NOTremoved+=seq[index]
            index+=1
            
    print("number of removed letters = ",index-len(NOTremoved))

    print(NOTremoved)
    return NOTremoved
#********************************************************

def make_words():
    query_seq = ""
    words_seq = []
    query_file = open("main seq.txt", "r")
    
    for line in query_file:
        query_seq= line
        
    
    print("number of letters = ",len(query_seq))
    print(query_seq)
    

    query=remove_repeats(query_seq)
    for j in range(len(query)):
        if(j+3 <= len(query)):
            words_seq.append(query[j: j+3])
        else:
            break
            
    print("number of words = ",len(words_seq))
    print(words_seq)
    query_file.close()
    return query_seq ,words_seq
#**********************************************************   
def neighbourhood():
    sequence , words = make_words()
    
    neighbour_words=[]
    s=""
    threshold = 13
    hspneighbours = []
    Protein_letters=['A','R','N','D','C','Q','E','G','H','I','L','K','M','F','P','S','T','W','Y','V','B','Z','X']

    for i in range(len(words)):
        each_word= words[i]
        for j in range(len(each_word)):
            for a in range(len(Protein_letters)):
                strr = list(each_word) #cut word into letters
                strr[j] = Protein_letters[a]
                for x in strr:  #convert list to string
                    s+=x
                if s not in neighbour_words:
                    neighbour_words.append(s)
                
                   #hsp
                    score = 0
                    for p in range(0,len(s)):
                       score+=int(matrix(s[p],each_word[p]))
                    
                    if score >= threshold:
                        hspneighbours.append((s , i, score))
                s=""
            #بمسك الكلمة وكل حرف فيها ابدله ب حروف بتوع البروتين       
      
    
    print("number of neighbour words = ",len(neighbour_words))
    print("number of high score neighbour words = ",len(hspneighbours))

    #for v in neighbour_words:
        #print (v)
    return neighbour_words , hspneighbours  #seed , query location, score  

#*************************************************************
def Hit(seed):
    hits = []
    f = open("seq1.txt", "r")
    db = f.read()
#    neighbours, seeds = neighbourhood()
    for k in range(len(seeds)):
        seed= seeds[k][0]
        for i in seed:
            for j in range(len(db)):
                if j + 3 < len(db):
                    if i == db[j] and db[j:j +3] == seed: #if the seed = word in database
                        hits.append((seed ,j))
                else:
                    break
                
    f.close()
    
    return hits    #seed , db location



#**************************************************************
neighbours, seeds = neighbourhood()
Hit(seeds)
