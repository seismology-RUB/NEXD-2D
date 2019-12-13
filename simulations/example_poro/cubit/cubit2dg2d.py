import cubit
import numpy
from math import sqrt

class mesh():
    def __init__(self):
        print('pass')
        # set main variables
        self.mesh_name='mesh'
        self.coord_name='coord'
        self.tri_name='tri'
        self.face_name='face'
        self.neighbor_name='neighbor'
        self.material_name='mat'
        self.matprop_name='matprop'
        self.absorb_name='absorb'
        self.free_name='free'
        self.tri='TRISHELL3'
        self.node='SPHERE'
        self.poro=False

        self.get_mesh()

        self.get_block()
        self.get_face()

        self.coord_write(self.coord_name)
        self.elem_write(self.mesh_name)
        self.mat_write(self.material_name)
        self.matprop_write(self.matprop_name)
        self.absorb_write(self.absorb_name)
        self.free_write(self.free_name)


    def get_block(self):
        ''' extract block information '''
        block_flag=[]
        block_mat=[]
        block_pml=[]
        block_pmlflag=[]
        ns_free=[]
        ns_absorb=[]
        material=[]
        # number of blocks
        blocks=cubit.get_block_id_list()

        # go through blocks
        for block in blocks:
            name=cubit.get_exodus_entity_name('block',block)
            type=cubit.get_block_element_type(block)
            print block,name,blocks,type

            if type == self.tri:
                # check if elastic
                if name.find('elastic') >=0:
                    print('found elastic block')
                    imat=1
                
                    # get number of attributes
                    values=name.split(" ")
                    print values
                    flag=int(values[1])
                    vp=float(values[2])
                    vs=float(values[3])
                    rho=float(values[4])
                    qp=float(values[5])
                    qs=float(values[6])
                    pml=float(values[7])

                    block_flag.append(flag)
                    block_pmlflag.append(pml)
                    block_mat.append(block)
                    material.append([imat,vp,vs,rho,qp,qs])
                elif name.find('poro') >=0:
                    print('found poroelastic block')
                    self.poro=True
                    imat=2
                    
                    # get number of attributes
                    values=name.split(" ")
                    print values
                    flag=int(values[1])
                    typeid=float(values[2])
                    rhos=float(values[3])
                    lambdau=float(values[4])
                    my=float(values[5])
                    phi=float(values[6])
                    kappa=float(values[7])
                    b=float(values[8])
                    invT=float(values[9])
                    invN=float(values[10])
                    rho1=float(values[11])
                    S1=float(values[12])
                    K1=float(values[13])
                    ny1=float(values[14])
                    rho2=float(values[15])
                    S2=float(values[16])
                    K2=float(values[17])
                    ny2=float(values[18])
                    fitting_n=float(values[19])
                    fitting_chi=float(values[20])
                    Sr1=float(values[21])
                    Sr2=float(values[22])

                    block_flag.append(flag)
                    block_mat.append(block)
                    material.append([imat,typeid,rhos,lambdau,my,phi,kappa,b,invT,invN,rho1,S1,K1,ny1,rho2,S2,K2,ny2,fitting_n,fitting_chi,Sr1,Sr2])
                elif name.find('pml') >=0:
                    ipml=1
                    values=name.split(" ")
                    flag=int(values[1])
                    block_pmlflag.append(flag)
                    block_pml.append(block)
                else:
                    print('error, no correct material in block',block)

        ns=cubit.get_nodeset_id_list()
        for n in ns:
            name=cubit.get_exodus_entity_name('nodeset',n)
            if name.find('free') >=0:
                print('found free surface nodes')
                ns_free.append(n)
            elif name.find('absorb') >=0:
                print('found absorb surface nodes')
                ns_absorb.append(n)
            elif name.find('vtk') >=0:
                pass
            else:
                print('error in boundaries',n)
            
        print('BLOCKFLAG',block_flag)
        self.block_flag=block_flag
        self.block_mat=block_mat
        self.block_pmlflag=block_pmlflag
        self.mat=material
        self.ns_free=ns_free
        self.ns_absorb=ns_absorb


    def get_mesh(self):
        ''' get tri mesh from cubit in format for dg2d '''
        print('Start extracting mesh for dg2d')

    def coord_write(self,coord_name):
        ''' write nodes file '''
        coord=open(coord_name,'w')
        print('writing '+coord_name)
        node_list=cubit.parse_cubit_list('node','all')
        ncoord=len(node_list)
        print('number of nodes:',str(ncoord))
        # write ncoord
        coord.write('%10i\n' % ncoord)
        #write coords
        for node in node_list:
            x,y,z=cubit.get_nodal_coordinates(node)
            txt=('%10i %20f %20f\n') % (node,x,y)
            coord.write(txt)
        coord.close()
        print('Ok')

    def mat_write(self,mat_name):
        nelem=cubit.get_tri_count()
        element=[ [0,0,0] for i in range(nelem)]
        for block,flag,pmlflag in zip(self.block_mat,self.block_flag,self.block_pmlflag):
            tris=cubit.get_block_tris(block)
            for tri in tris:
                element[tri-1] = [tri,flag,pmlflag]
        mat=open(mat_name,'w')
        print('Writing '+mat_name+'.....')
        for i in range(nelem):
            mat.write(('%10i %10i %10i\n') % (element[i][0],element[i][1],element[i][2]))
        mat.close()
        print('Ok')

    def elem_write(self,mesh_name):
        meshfile=open(mesh_name,'w')
        print('Writing '+mesh_name+'.....')
        nelem=cubit.get_tri_count()
        print('number of elements:',str(nelem))
        meshfile.write(str(nelem)+'\n')
        num_write=0
        temp_tri=[]
        trilength=1
        for block,flag in zip(self.block_mat,self.block_flag):
            trilength += len(cubit.get_block_tris(block))
        tri_vp=range(trilength)
        tri_vs=range(trilength)
        tri_rho=range(trilength)
        tri_qp=range(trilength)
        tri_qs=range(trilength)
        tri_block=range(trilength)
        for block,flag in zip(self.block_mat,self.block_flag):
            tris=cubit.get_block_tris(block)
            name=cubit.get_exodus_entity_name('block',block)
            type=cubit.get_block_element_type(block)
            for tri in tris:
                nodes=cubit.get_connectivity('Tri',tri)
                temp_tri.append(tri)
                values=name.split(" ")
                tri_vp[tri]=float(values[2])
                tri_vs[tri]=float(values[3])
                tri_rho[tri]=float(values[4])
                tri_qp[tri]=float(values[5])
                tri_qs[tri]=float(values[6])
        temp_tri.sort()
        for tri in temp_tri:
            nodes=cubit.get_connectivity('Tri',tri)
            txt=('%10i ')% tri
            txt=txt+('%10i %10i %10i\n')% nodes[:]
            # This is later needed for the inversion feature:
            #txt=txt+('%10i %10i %10i')% nodes[:]
            #txt=txt+(' %9.1f %9.1f %9.1f %5i %5i\n') % (tri_vp[tri],tri_vs[tri],tri_rho[tri],tri_qp[tri],tri_qs[tri])
            meshfile.write(txt)


        meshfile.close()
        print('Ok')

    def matprop_write(self,matprop_name):
        if self.poro:
            for material in self.mat:
                if len(material)<22:
                    raise ValueError('Please only define poroelastic materials and do not include other materials.')

        matpropfile=open(matprop_name,'w')
        print('writing '+matprop_name)
        nmat=len(self.block_mat)
        print('number of materials: ',str(nmat))
        matpropfile.write(str(nmat)+'\n')
        for i,block in enumerate(self.block_mat):
            txt=str(block)
            if self.poro:
                txt= (txt+' '+str(int(self.mat[i][1]))+' '+str(self.mat[i][2])+' '+str(self.mat[i][3])+' '+str(self.mat[i][4])+' '
                      +str(self.mat[i][5])+' '+str(self.mat[i][6])+' '+str(self.mat[i][7])+' '+str(self.mat[i][8])+' '+str(self.mat[i][9])+' '
                      +str(self.mat[i][10])+' '+str(self.mat[i][11])+' '+str(self.mat[i][12])+' '+str(self.mat[i][13])+' '+str(self.mat[i][14])+' '
                      +str(self.mat[i][15])+' '+str(self.mat[i][16])+' '+str(self.mat[i][17])+' '+str(self.mat[i][18])+' '+str(self.mat[i][19])+' '
                      +str(self.mat[i][20])+' '+str(self.mat[i][21]))
            else:
                txt= (txt+' '+str(self.mat[i][0])+' '+str(self.mat[i][1])+' '+str(self.mat[i][2])+' '+str(self.mat[i][3])+' '
                      +str(self.mat[i][4])+' '+str(self.mat[i][5]))
            matpropfile.write(txt+'\n')
        
    def get_face(self):
        for n in self.ns_free:
            nodes_free=cubit.get_nodeset_nodes(n)
            print(nodes_free)

    def absorb_write(self,absorb_name):
        absorbfile=open(absorb_name,'w')
        print('writing '+absorb_name)

        if self.ns_absorb!=[]:
            for n in self.ns_absorb:
                nodes_absorb=cubit.get_nodeset_nodes(n)
            nabs=len(nodes_absorb)
            absorbfile.write(str(nabs)+'\n')
            print('Number of absorbing nodes ',nabs)
            for node in nodes_absorb:
                absorbfile.write(str(node)+'\n')
        else:
            absorbfile.write(str(0)+'\n')
            print('Number of absorbing nodes ',0)

    def free_write(self,free_name):
        freefile=open(free_name,'w')
        print('writing '+free_name)

        if self.ns_free!=[]:
            for n in self.ns_free:
                nodes_free=cubit.get_nodeset_nodes(n)
            nfree=len(nodes_free)
            freefile.write(str(nfree)+'\n')
            print('Number of free nodes ',nfree)
            for node in nodes_free:
                freefile.write(str(node)+'\n')
        else:
            freefile.write(str(0)+'\n')
            print('Number of free nodes ',0)


if __name__ == '__main__':
    mesh()
