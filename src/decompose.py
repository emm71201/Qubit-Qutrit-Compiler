import numpy
from node import Node
import scipy
from utils import *
import tqdm

# log z = log |z| + i Arg z

def bfs(root):

    queue = [root]
    visited = []

    while len(queue) > 0:
        curr = queue.pop(0)
        if curr is not None and curr not in visited:
            #print(curr, len(queue), len(visited))

            if curr.left is not None:
                queue.append(curr.left)
            if curr.right is not None:
                queue.append(curr.right)

            visited.append(curr)

    #print("#"*70)
    #print(curr, len(queue), len(visited))




# decomposition is a breadth first traversal
def decompose(unitary, register):
    regmap = ["qubit" if item == 2 else "qutrit" for item in register]
    root = Node(unitary, 0, None)

    queue = [root]
    visited = []


    pbar = tqdm.tqdm(total=len(queue), desc="Decomposing the Unitary Matrix")
    while len(queue) > 0:
        last_len_queue = len(queue)
        curr = queue.pop(0)
        if curr is not None and curr not in visited:
            # do stuff

            #print(f"Len queue = {len(queue)}")
            if register[-1] == 2 and register[-2] == 2:
                if curr.type is None or curr.type == -1 and curr.dim > 4:
                    curr.decompose(regmap[curr.level])
            else:
                if curr.type is None or curr.type == -1 and curr.dim != 2:
                    curr.decompose(regmap[curr.level])
                    print("lalala")

            for key, child in curr.children.items():
                if child is not None and not child in visited:
                    if child.is_identity():
                        child.parent.children[key] = None
                        print("Identity found")
                    else:
                        queue.append(child)

            visited.append(curr)

        # pbar.update(1)
        # pbar.total = len(queue)
        if len(queue) > last_len_queue:
            pbar.reset(0)
            pbar.total = len(queue)
            #pbar.postfix = f"Length of Decomposition: {len(queue)}"
            pbar.update(0)
        else:
            #pbar.postfix = f"Length of Decomposition: {len(queue)}"
            pbar.update(1)
    pbar.close()


    return root



if __name__=="__main__":
    matrix = numpy.genfromtxt("matrices/s108fft-108dim.csv", delimiter=",", dtype=numpy.complex_)
    register = [3,3,3,2,2]

    root = decompose(matrix, register)
    child = root.children[0].children[0].children[0]
    child2 = root.children[0].children[0].children[10]
    print(child.dim)
    print(child.op)
    print("#"*70)
    print(child2.op)
    # root.chi
    # root.decompose(type="qutrit", verify=True)
    # child1 = root.children[1]
    # print(root.children[22].op / numpy.pi)
    # child1.decompose("qutrit")

    # a = Node("A")
    # b = Node("B")
    # c = Node("C")
    # a.left = b
    # a.right = c
    # d = Node("D")
    # e = Node("E")
    # b.left = d
    # b.right = e
    # f = Node("F")
    # g = Node("G")
    # e.left = f
    # e.right = g
    #
    # bfs(a)


