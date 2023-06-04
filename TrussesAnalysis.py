import numpy as np

class Point:#结点类
    def __init__(self, x, y,index,load = [0,0],pin = False):
        self.x = x
        self.y = y
        self.index = index
        self.load = load
        self.pin = pin
        
class Element:#杆件类
    def __init__(self, point1, point2,EA=1):
        self.point1 = point1
        self.point2 = point2
        self.EA = EA
    def length(self):
        return np.sqrt((self.point2.x - self.point1.x)**2 + (self.point2.y - self.point1.y)**2)
    def transform_Mat(self):
        theta = np.arctan2(self.point2.y - self.point1.y, self.point2.x - self.point1.x)
        transform_mat = np.array([[np.cos(theta), np.sin(theta),0,0],
                                 [0,0,np.cos(theta),np.sin(theta)]])
        return transform_mat
    def stiffness_Mat(self):
        return np.array([[1,-1],[-1,1]])*self.EA/self.length()
    def stiffness_Mat_transed(self):
        '''
        返回转换后的刚度矩阵T.T()x K x T
        '''

        local_stiffness_mat = np.array([[1,-1],[-1,1]])*self.EA/self.length()
        transform_mat = self.transform_Mat()
        return transform_mat.transpose()@local_stiffness_mat@transform_mat

    def location_vec(self):
        """
        返回定位向量

        :return: numpy array containing location vectors.
        """
        return np.array([self.point1.index*2,self.point1.index*2+1,self.point2.index*2,self.point2.index*2+1])
    
    
    
class Truss:
    def __init__(self,nodes:list,elements:list) -> None:
        self.nodes = nodes
        self.elements = elements
        self.K = np.zeros((len(self.nodes)*2,len(self.nodes)*2)) #初始化总刚度矩阵，维度为（结点数*2，结点数*2）
        

    def integrate_K(self):
        for e in self.elements:
            k_e = e.stiffness_Mat_transed()
            location_vec = e.location_vec()
            for row in range(4):
                for col in range(4):
                    if k_e[row,col] != 0:
                        self.K[location_vec[row],location_vec[col]] += k_e[row,col]
    def get_Load_vec(self):#P
        load_vec = np.zeros((len(self.nodes)*2,1))
        for node in self.nodes:
            if node.load!=[0,0]:
                load_vec[node.index*2] = node.load[0]
                load_vec[node.index*2+1] = node.load[1]
        return load_vec
    def restrict(self):#处理铰接点约束:置0置1法
        for node in self.nodes:
            if node.pin:
                self.K[node.index*2,:] = 0
                self.K[:,node.index*2] = 0
                self.K[node.index*2+1,:] = 0
                self.K[:,node.index*2+1] = 0
                self.K[node.index*2,node.index*2] = 1
                self.K[node.index*2+1,node.index*2+1] = 1
                
    def calc(self):
        self.integrate_K()
        self.restrict()
        Deltas = np.linalg.inv(self.K)@self.get_Load_vec()
        inner_forces = []
        for e in self.elements:
            location_vec = e.location_vec()
            delta_e = e.transform_Mat()@ Deltas[location_vec]
            force = e.stiffness_Mat()@delta_e
            inner_forces.append(force)
        return inner_forces
def test():
    nodes = [Point(0,0,0,pin=True,load=[0,0]),Point(4,0,1,pin=True),Point(8,0,2,pin=True),\
        Point(2,3,3,load=[0,-10]),Point(6,3,4,load=[0,-10])]
    elements = [Element(nodes[0],nodes[1]),Element(nodes[1],nodes[2]),Element(nodes[0],nodes[3]),\
        Element(nodes[3],nodes[1]),Element(nodes[4],nodes[1]),Element(nodes[4],nodes[2]),Element(nodes[3],nodes[4])]
    t = Truss(nodes=nodes,elements=elements)
    
    res = t.calc()
    for i in res:
        print(i[0])
test()