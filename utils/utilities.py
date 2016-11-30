import numpy as np

class Utilities:
    def str5sf(self,number):
        num_str = str(number)
        while len(num_str) <= 5:
            num_str = num_str+'0'
        return num_str

    def my_reshape(self,M):
        # change order to nc,mq,ms,mu,nboot
        dim = np.shape(M)
        M2 = np.zeros([dim[0],dim[3],dim[2],dim[1],dim[4]])
        for i in range(dim[0]):
            for j in range(dim[1]):
                for k in range(dim[2]):
                    for l in range(dim[3]):
                        for m in range(dim[4]):
                            M2[i,l,k,j,m] = M[i,j,k,l,m]
        return M2

    def ZAZP(self,kin,basis):
        ZA = [0.71273,0.74404]
        """for i in range(10):
            if i==0:
            elif i==1:
                elif i==2 or 9:
                    elif i==3 or 8:
                        elif i==4 or 7:
                            elif i==5 or 6:"""
        ZP = [0.57732,0.57773] #stopgap
        
        return ZA[basis], ZP[basis]

    def matmul(self,A,B):
        n = A.shape[0]
        m = B.shape[1]
        l1 = A.shape[1]
        l2 = B.shape[0]

        if l1==l2:
            C = np.zeros([n,m])
            for i in range(n):
                for j in range(m):
                    sum = 0
                    for k in range(l1):
                        sum = sum + A[i][k]*B[k][j]
                    C[i,j] = sum

            return C
        else:
            print "Error, matrices dimensions are wrong"
            print "Error, matrices dimensions are wrong"
            return -1

    def transpose(self,A):
        l1 = A.shape[0]
        l2 = A.shape[1]
        
        B=np.zeros([l2,l1])
        for i in range(l1):
            for j in range(l2):
                B[j,i] = A[i,j]
        
        return B

    def diagonalize(self,A):
        dim = A.shape[0]
        if(dim != A.shape[1]):
            print "Matrix not square"
            print "Matrix not square"
            return -1
        else:
            B = np.zeros([dim,1])
            for i in range(dim):
                B[i] = A[i,i]
            return B


