import numpy as np



class Node:#结点类
    def __init__(self, x, y,index,load = [0,0,0],support = 0):
        """
        load: 荷载（水平力，竖直力，弯矩（逆时针为正））
        support: 0: no support; 
        1 : pin
        2: horizantal roller
        3: vertical roller
        4: fixed
        """
        self.x = x
        self.y = y
        self.index = index
        self.load = load
        self.support = support
class Element:#杆件类
    def __init__(self, node1, node2,E=1,A=1,I=1):
        self.node1 = node1
        self.node2 = node2
        self.E = E
        self.A = A
        self.I = I
    def length(self):
        return np.sqrt((self.node2.x - self.node1.x)**2 + (self.node2.y - self.node1.y)**2)
    def transform_Mat(self):
        theta = np.arctan2(self.point2.y - self.point1.y, self.point2.x - self.point1.x)
        transfrom_mat = np.array([[np.cos(theta), np.sin(theta),0,0,0,0],],
                                 [-np.sin(theta), np.cos(theta),0,0,0,0],
                                 [0,0,0,np.cos(theta),np.sin(theta),0],
                                 [0,0,0,-np.sin(theta),np.cos(theta),0],
                                 [0,0,0,0,0,1])
        return transfrom_mat
    def stiffness_Mat(self):
        l = self.length()
        E = self.E
        A = self.A
        I = self.I
        return np.array(
            [
                [E*A/l,0,0,-E*A/l,0,0],
                [0,12*E*I/l**3,6*E*I/l**2,0,-12*E*I/l**3,6*E*I/l**2],
                [0,6*E*I/l**2,4*E*I/l,0,-6*E*I/l**2,2*E*I/l],
                [-E*A/l,0,0,E*A/l,0,0],
                [0,-12*E*I/l**3,-6*E*I/l**2,0,12*E*I/l**3,-6*E*I/l**2],
                [0,6*E*I/l**2,2*E*I/l,0,-6*E*I/l**2,4*E*I/l]
            ]
        )
    def stiffness_Mat_transed(self):
        local_stiffness_mat = self.stiffness_Mat()#局部刚度矩阵
        transform_mat = self.transform_Mat()
        return transform_mat.transpose()@local_stiffness_mat@transform_mat
    
    def location_vec(self):
        """
        返回定位向量(1*6)
        """
        return np.array([
            self.node1.index*3,self.node1.index*3+1,self.node1.index*3+2,
            self.node2.index*3,self.node2.index*3+1,self.node2.index*3+2
        ])


class Frame:
    def __init__(self,nodes:list,elements:list) -> None:
        self.nodes = nodes
        self.elements = Element
        self.stiffnessMat = np.zeros((len(nodes)*3,len(nodes)*3))#初始化总刚度矩阵
    def integrate_K(self):#组装总刚度矩阵
        for e in self.elements:
            k_e = e.stiffness_Mat_transed()
            location_vec = e.location_vec()
            for row in range(6):
                for col in range(6):
                    if k_e[row,col] !=0:
                        self.stiffnessMat[location_vec[row],location_vec[col]] += k_e[row,col]
    def get_Load_vec(self):#荷载向量
        load_vec = np.zeros((len(self.nodes)*3,1))
        for node in self.nodes:
            if node.load!=[0,0,0]:
                load_vec[node.index*3] = node.load[0]
                load_vec[node.index*3+1] = node.load[1]
                load_vec[node.index*3+2] = node.load[2]
        return load_vec