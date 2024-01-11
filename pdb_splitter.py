import os
import sys

#split our file, call with 0 for struc, 1 for cov 

ter_state = 0 

def check_input(args):
    """Checks whether to read from stdin/file and validates user input/options.
    """

    # Defaults
    option = ''
    fh = sys.stdin  # file handle
    option = args[0][1]
    fh = open(args[1], 'r')
    return (fh, option )

def run(fhandle, option):

    #ignore HET and other line starts: 
    records = ('ATOM', 'ANISOU', 'TER')
    for line in fhandle:
        if line.startswith(records):
            
            if (option == 0) and (ter_state==0)  : 
                yield line
                
                if line.startswith('Ter') :
                    break #break once we get to first Ter as option 0
            
            #need to check for Ter after store line starting with Ter due to structure 
            #of pdb files, TER line belongs to structure. 
            if line.startswith('Ter')and (ter_state==0) :
                ter_state =1
                continue 
            
            if (option == 1 ) and (ter_state==1): 
                yield line
            else:
                continue 
            
delete_residue_by_name = run

def main():
    # Check Input TODO fix input arg selection 
    pdbfh, resname_set , idx= check_input(sys.argv[1:])

    # Do the job
    new_pdb = run(pdbfh, resname_set, idx)

    try:
        _buffer = []
        _buffer_size = 5000  # write N lines at a time
        for lineno, line in enumerate(new_pdb):
            if not (lineno % _buffer_size):
                sys.stdout.write(''.join(_buffer))
                _buffer = []
            _buffer.append(line)

        sys.stdout.write(''.join(_buffer))
        sys.stdout.flush()
    except IOError:
        # This is here to catch Broken Pipes
        # for example to use 'head' or 'tail' without
        # the error message showing up
        pass

    # last line of the script
    # We can close it even if it is sys.stdin
    pdbfh.close()
    sys.exit(0)


if __name__ == '__main__':
    main()