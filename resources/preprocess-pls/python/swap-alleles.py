import sys

def main():
    # Read directly from standard input (the pipe)
    for line in sys.stdin:
        if line.startswith('#'):
            sys.stdout.write(line)
            continue
            
        cols = line.rstrip('\n').split('\t')
        info = cols[7]
        
        # Parse out the target annotations
        bpos, bmap, bref, balt = None, None, None, None
        
        for item in info.split(';'):
            if item.startswith('BPOS='): bpos = item[5:]
            elif item.startswith('BMAP='): bmap = item[5:]
            elif item.startswith('BREF='): bref = item[5:]
            elif item.startswith('BALT='): balt = item[5:]
                
        # Perform the column swap if the tags were found
        if bpos is not None: cols[1] = bpos
        if bmap is not None: cols[2] = bmap
        if bref is not None: cols[3] = bref
        if balt is not None: cols[4] = balt
            
        sys.stdout.write('\t'.join(cols) + '\n')

if __name__ == "__main__":
    main()
