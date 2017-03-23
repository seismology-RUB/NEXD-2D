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

                    block_flag.append(flag)
                    block_mat.append(block)
                    material.append([imat,vp,vs,rho,qp,qs])
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
            else:
                print('error in boundaries',n)
            
        print('BLOCKFLAG',block_flag)
        self.block_flag=block_flag
        self.block_mat=block_mat
        self.block_pml=block_pml
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
        for block,flag in zip(self.block_mat,self.block_flag):
            tris=cubit.get_block_tris(block)
            for tri in tris:
                element[tri-1] = [tri,flag,0]

        for block,flag in zip(self.block_pml,self.block_pmlflag):
            tris=cubit.get_block_tris(block)
            for tri in tris:
                element[tri-1][2] = flag

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
        for block,flag in zip(self.block_mat,self.block_flag):
            tris=cubit.get_block_tris(block)
            for tri in tris:
                nodes=cubit.get_connectivity('Tri',tri)
                temp_tri.append(tri)
#                txt=('%10i ')% tri
#                txt=txt+('%10i %10i %10i\n')% nodes[:]
 #               meshfile.write(txt)

        temp_tri.sort()
        for tri in temp_tri:
            nodes=cubit.get_connectivity('Tri',tri)
            txt=('%10i ')% tri
            txt=txt+('%10i %10i %10i\n')% nodes[:]
            meshfile.write(txt)


        meshfile.close()
        print('Ok')

    def matprop_write(self,matprop_name):
        matpropfile=open(matprop_name,'w')
        print('writing '+matprop_name)
        nmat=len(self.block_mat)
        print('number of materials: ',str(nmat))
        matpropfile.write(str(nmat)+'\n')
        for i,block in enumerate(self.block_mat):
#            print(block)
#            print(block,self.mat[i][:])
            txt=str(block)
#            txt=txt+' '+str(self.mat[i][0])+' '+str(self.mat[i][1])+' '+str(self.mat[i][2])+' '+str(self.mat[i][3])
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
        
        for n in self.ns_absorb:
            nodes_absorb=cubit.get_nodeset_nodes(n)
        nabs=len(nodes_absorb)
        absorbfile.write(str(nabs)+'\n')
        print('Number of absorbing nodes ',nabs)

        for node in nodes_absorb:
            absorbfile.write(str(node)+'\n')

    def free_write(self,free_name):
        freefile=open(free_name,'w')
        print('writing '+free_name)
        
        for n in self.ns_free:
            nodes_free=cubit.get_nodeset_nodes(n)
        nfree=len(nodes_free)
        freefile.write(str(nfree)+'\n')
        print('Number of free nodes ',nfree)

        for node in nodes_free:
            freefile.write(str(node)+'\n')

        


if __name__ == '__main__':
    mesh()
