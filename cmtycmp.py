import sys


import pcd.cmty as cmty
import pcd.cmtycmp
import pcd.tcmty as tcmty


def dyncmp(c1a, c1b, c2a, c2b, cmpfunc=pcd.cmtycmp.nmi):
    nc1t1 = c1t1.nodecmtys_onetoone()
    nc1t2 = c1t2.nodecmtys_onetoone()
    nc2t1 = c2t1.nodecmtys_onetoone()
    nc2t2 = c2t2.nodecmtys_onetoone()

    c1_transitions = { }
    c2_transitions = { }
    for n in nc1t1:
        c1_transitions[n] = (nc1t1[n], nc1t2[n])
        c2_transitions[n] = (nc2t1[n], nc2t2[n])
        #print c1_transitions[n], c2_transitions[n]

    transition_cmtys1 = cmty.Communities.from_nodecmtys(c1_transitions)
    transition_cmtys2 = cmty.Communities.from_nodecmtys(c2_transitions)

    #print transition_cmtys1.cmtysizes(), transition_cmtys2.cmtysizes()
    #print pcd.cmtycmp.confusion(transition_cmtys1, transition_cmtys2)
    value = cmpfunc(transition_cmtys1, transition_cmtys2)
    return value

def dyncmp_series(tcmtys1, tcmtys2, dt, cmpfunc):

    print "#t1 t2 value"
    ts = list(tcmtys1)
    if   dt  > 0:             endpoint = -dt
    elif dt == 0:    endpoint = None
    for i, t1 in enumerate(ts[:endpoint]):
        t2 = ts[i+dt]
        c1t1 = tcmtys1[t1]
        c1t2 = tcmtys1[t2]
        c2t1 = tcmtys2[t1]
        c2t2 = tcmtys2[t2]

        value = dyncmp(c1t1, c1t2, c2t1, c2t2, cmpfunc=cmpfunc)

        print t1, t2, value


if __name__ == '__main__':
    import argparse

    parser = argparse.ArgumentParser(description='Process some integers.')
    parser.add_argument('input1', help='Input file #1')
    parser.add_argument('input2', help='Input file #2')
    parser.add_argument('--format', default='auto',
                        help='Format for input files: auto, '
                             'tcommlist, tmatrix (default: %(default)s).')
    parser.add_argument('--dt', default=1, type=int,
                        help='time window (default=%(default)s)')
    parser.add_argument('--cmp', default='nmi',
                        help='comparison function: vi, nmi, '
                             'mutual_information, jaccard, F1, '
                             'etc (default: %(default)s)')

    args = parser.parse_args()


    reader = getattr(tcmty.TemporalCommunities, 'from_'+args.format)
    tcmtys1 = reader(sys.argv[1])
    tcmtys2 = reader(sys.argv[2])

    dyncmp(tcmtys1, tcmtys2, dt=args.dt,
           cmpfunc=getattr(pcd.cmtycmp, args.cmp))
