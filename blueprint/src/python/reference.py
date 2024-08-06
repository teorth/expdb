# Code to read and interpret bibtex references 

class Reference:
    def __init__(self, label, reftype, fields):
        self.label = label
        self.reftype = reftype
        self.entries = fields
        
    def __repr__(self):
        title = self.entries["title"] if "title" in self.entries else "[No title]"
        return f'{self.author()} ({self.year()}) {title}'
    
    def author(self):
        return self.entries['author'] if 'author' in self.entries else 'Unknown author'
    
    def year(self):
        return self.entries['year'] if 'year' in self.entries else 'Unknown date'
    
    
    def is_classical(self):
        return self.label == 'Classical'
    
    def is_derived(self):
        return self.label == 'Derived'
    
    def is_conjectured(self):
        return self.label == 'Conjectured'
    
    def is_trivial(self):
        return self.label == 'Trivial'
        
    def is_literature(self):
        return not (self.is_classical() or self.is_derived() or self.is_conjectured())
        
    # Static functions --------------------------------------------------------
    # Constructor for a (author, year) reference
    def make(author, year):
        label = f'temp_{author}_{year}'
        reftype = 'article'
        fields = {'author': author, 'year': year}
        return Reference(label, reftype, fields)    
    
    # Constructor for a classical reference
    def classical():    
        return Reference('Classical', 'Classical', { 'year': -1 })
    
    # Constructor for a derived result
    def derived(year):
        return Reference('Derived', 'Derived', { 'year': year })
    
    # Constructor for a conjectured result
    def conjectured():
        return Reference('Conjectured', 'Conjectured', { 'year': float('inf') })
    
    def trivial():
        return Reference('Trivial', 'Trivial', { 'year': -1 })
    
    # Returns the maximum year among a list of references. This is used to keep 
    # track of when a lemma could have been proved. 
    def max_year(refs):
        if len(refs) == 0: return -1
        if any(r.year() == 'Unknown date' for r in refs): return 'Unknown date'
        return max(r.year() for r in refs)
    
class Reference_Manager:
    
    def __init__(self, path):
        self.path = path
        self.refs = {}
    
    # Read the bibtex file and store references as Reference objects. Assumes 
    # label of the references are unique. 
    def load(self):
        file = open(self.path, 'r')
        lines = file.readlines()
        file.close()
        
        # parse file 
        s = ' '.join(l.strip() for l in lines)
        items = s.split('@')
        
        for item in items:
            if len(item) == 0: continue
            ref = self._parse(item)
            if len(ref.label) > 0:
                self.refs[ref.label] = ref
        
    def find_all(self, author='Any', year='Any'):
        refs = []
        for k in self.refs:
            ref = self.refs[k]
            if (author == 'Any' or author in ref.entries["author"]) and \
                (year == 'Any' or ref.entries["year"] == year):
                    refs.append(ref)
        return refs
    
    # Get a reference
    def get(self, label):
        if label in self.refs:
            return self.refs[label]
        raise ValueError(f'Cannot find {label}')
    
    #--------------------------------------------------------------------------
    # Internal methods 
    
    # Given a string of the form 'A{B}C{{D}}', returns B{D}. 
    def _extract(self, item, open_brace, close_brace):
        chs = []
        level = 0
        for c in item:
            if c == close_brace: level -= 1
            if level >= 1: chs.append(c)
            if c == open_brace: level += 1
        return ''.join(chs)
    
    # Split by delimiter, but only when it is nested 'level' deep
    def _split_level(self, item, open_brace, close_brace, delim, level):
        split = [] # contains all words
        chs = [] # current word
        l = 0
        for c in item:
            if c == close_brace: l -= 1
            if l == level and c == delim:
                if len(chs) > 0:
                    split.append(''.join(chs))
                    chs.clear()
            else:
                chs.append(c)
            if c == open_brace: l += 1
        
        if len(chs) > 0:
            split.append(''.join(chs))
        return split
    
    def _parse_year(self, s):
        return int(s.split('/')[0]) # year is (hopefully) stored as an integer
        
    def _parse_authors(self, s):
        authors = [p.strip() for p in s.split('and')]
        for i in range(len(authors)):
            a = authors[i]
            if ',' in a:
                authors[i] = a.split(',')[0]
            else:
                authors[i] = a.split(' ')[-1]
        return '--'.join(authors)
        
    # Parses a string representing a reference item, returns as a Reference object
    def _parse(self, item):
        
        typ = item.split('{')[0].lower()
        ref = self._extract(item, '{', '}')
        entries = self._split_level(ref, '{', '}', ',', 0)
        label = entries[0].strip()
        
        fields = {}
        for i in range(1, len(entries)):
            e = entries[i].strip()
            if len(e) == 0: continue
            kvp = self._split_level(e, '{', '}', '=', 0)
            if len(kvp) >= 2:
                key = kvp[0].strip().lower()
                value = self._extract(kvp[1].strip(), '{', '}')
                if len(key) == 0 or len(value) == 0: continue
                if key == 'year': value = self._parse_year(value)
                elif key == 'author': value = self._parse_authors(value)
                fields[key] = value
                
        return Reference(label, typ, fields)
        
        