class Colors:
    RED = '\33[31m'
    GREEN = '\33[32m'
    BLUE = '\33[34m'
    YELLOW = '\33[33m'
    END = '\033[0m'


# -------------------------------------------------------------------------------------


class Atom:
    def __init__(self, atom_name, element_symbol, x_coordinate, y_coordinate, z_coordinate):
        try:
            self.atom_name = str(atom_name).strip()
            self.element_symbol = str(element_symbol).strip()
            self.x_coordinate = float(x_coordinate)
            self.y_coordinate = float(y_coordinate)
            self.z_coordinate = float(z_coordinate)
        except:
            raise ValueError


# -------------------------------------------------------------------------------------


class Termination_symbol:
    def __init__(self, residue_name, residue_sequence_number):
        try:
            self.residue_name = str(residue_name).strip()
            self.residue_sequence_number = int(residue_sequence_number)
        except:
            raise ValueError


# -------------------------------------------------------------------------------------


class Residue:
    def __init__(self, residue_name, residue_sequence_number):
        try:
            self.residue_name = str(residue_name).strip()
            self.residue_sequence_number = int(residue_sequence_number)
            self.ATOMS = []
        except:
            raise ValueError

    def push_atom(self, next_atom):
        self.ATOMS.append(next_atom)


# -------------------------------------------------------------------------------------


class Second_order_structure:
    def __init__(self, chain_identifier, structure_identifier, residues_list):
        try:
            self.chain_identifier = str(chain_identifier).strip()
            self.structure_identifier = str(structure_identifier).strip()
            self.RESIDUES = residues_list
        except:
            raise ValueError


# -------------------------------------------------------------------------------------


class Helix(Second_order_structure):
    def __repr__(self):
        residues_names = ', '.join([residue.residue_name for residue in self.RESIDUES])
        return f'{Colors.GREEN}CHAIN: {self.chain_identifier},{Colors.BLUE} HELIX: {self.structure_identifier} --> {Colors.END}{residues_names}'


# -------------------------------------------------------------------------------------


class Sheet(Second_order_structure):
    def __repr__(self):
        residues_names = ', '.join([residue.residue_name for residue in self.RESIDUES])
        return f'{Colors.GREEN}CHAIN: {self.chain_identifier},{Colors.YELLOW} SHEET: {self.structure_identifier} --> {Colors.END}{residues_names}'


# -------------------------------------------------------------------------------------


class Chain:
    def __init__(self, chain_name):
        try:
            self.chain_name = str(chain_name).strip()
            self.RESIDUES = []
            self.HELICES = []
            self.SHEETS = []
            self.helices_indexes = []
            self.sheets_indexes = []
        except:
            raise ValueError

    def __repr__(self):
        chain_arrow = f'{Colors.GREEN}CHAIN: {self.chain_name} --> {Colors.END}'
        residues_names = ''

        for i in range(len(self.RESIDUES)):
            residues_names += self.RESIDUES[i].residue_name
            if ((i + 1) % 10 == 0):
                residues_names += '\n'
                residues_names += ''.ljust(len(chain_arrow) - 9)
            else:
                residues_names += ' '

        return f'{chain_arrow}{residues_names}'

    def push_residue(self, next_residue):
        if type(next_residue) == Residue:
            self.RESIDUES.append(next_residue)
        else:
            raise ValueError

    def push_helix_indexes(self, helix_identifier, start_index, end_index):
        self.helices_indexes.append( (helix_identifier, start_index, end_index) )

    def push_helix(self, chain_identifier, helix_identifier, residues_list):
        self.HELICES.append(Helix(chain_identifier, helix_identifier, residues_list))

    def push_sheet_indexes(self, sheet_identifier, start_index, end_index):
        self.helices_indexes.append( (sheet_identifier, start_index, end_index) )

    def push_sheet(self, chain_identifier, sheet_identifier, residues_list):
        self.SHEETS.append(Sheet(chain_identifier, sheet_identifier, residues_list))

    def get_residues_names(self):
        return [residue.residue_name for residue in self.RESIDUES]

    def get_residues_number(self):
        return len(self.RESIDUES)

    def get_atoms_number(self):
        atoms_number = 0
        for residue in self.RESIDUES:
            atoms_number += len(residue.ATOMS)

        return atoms_number


# -------------------------------------------------------------------------------------


class Structure:
    def __init__(self, pdb_file_path):
        self.data_correctnes_flag = True
        self.structure_name = ""
        self.CHAINS = []

        try:
            pdb_file = open(pdb_file_path, "r")
            helices_indexes_raw = []
            sheets_indexes_raw = []

            prev_chain_name = ''
            prev_residue_number = -1
            pdb_file_line_counter = 0

# -------------------------------------------------------------------------------------
# Główna pętla parsera
# Oznaczenia według: https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/tutorials/pdbintro.html
# -------------------------------------------------------------------------------------

            try:
                for DATA in pdb_file:
                    pdb_file_line_counter += 1
                    
                    if DATA[0 : 5] == 'HELIX':         # HELIX_IDENT                 CHAIN_IDENT                 RES_SEQ_NUM         RES_SEQ_NUM
                        helices_indexes_raw.append( (str(DATA[11 : 14]).strip(), str(DATA[19 : 20]).strip(), int(DATA[21 : 25]), int(DATA[33 : 37])) )

                    elif DATA[0 : 5] == 'SHEET':      # SHEET_IDENT                 CHAIN_IDENT                 RES_SEQ_NUM         RES_SEQ_NUM
                        sheets_indexes_raw.append( (str(DATA[11 : 14]).strip(), str(DATA[21 : 22]).strip(), int(DATA[22 : 26]), int(DATA[33 : 37])) )

                    elif DATA[0 : 4] == 'ATOM' or DATA[0 : 6] == 'HETATM':
                                       # ATOM_NAME      ELEM_SYMBOL    X_COORD        Y_COORD        Z_COORD
                        next_atom = Atom(DATA[12 : 16], DATA[76 : 78], DATA[30 : 38], DATA[38 : 46], DATA[46 : 54])

                        if prev_chain_name == DATA[21 : 22]: # CHAIN_IDENT
                            if prev_residue_number == DATA[22 : 26]: # RES_SEQ_NUM
                                self.CHAINS[-1].RESIDUES[-1].push_atom(next_atom)
                            else:                    # RES_NAME       RES_SEQ_NUM
                                next_residue = Residue(DATA[17 : 20], DATA[22 : 26])
                                self.CHAINS[-1].push_residue(next_residue)
                                self.CHAINS[-1].RESIDUES[-1].push_atom(next_atom)
                                prev_residue_number = DATA[22 : 26] # RES_SEQ_NUM

                        else:
                            next_chain = Chain(DATA[21 : 22]) # CHAIN_IDENT
                                                 # RES_NAME       RES_SEQ_NUM
                            next_residue = Residue(DATA[17 : 20], DATA[22 : 26])
                            self.CHAINS.append(next_chain)
                            self.CHAINS[-1].push_residue(next_residue)
                            self.CHAINS[-1].RESIDUES[-1].push_atom(next_atom)
                            prev_chain_name = DATA[21 : 22] # CHAIN_IDENT
                            prev_residue_number = DATA[22 : 26] # RES_SEQ_NUM

                    elif DATA[0 : 3] == 'TER':                           # RES_NAME       RES_SEQ_NUM
                        self.CHAINS[-1].RESIDUES.append(Termination_symbol(DATA[17 : 20], DATA[22 : 26]))

                pdb_file.close()
                
            except:
                self.data_correctnes_flag = False
                print(Colors.RED + f'Plik \"{pdb_file_path}\" zawiera blad w {pdb_file_line_counter} linii :O' + Colors.END)

# -------------------------------------------------------------------------------------
# Usunięcie fragmentów pochopnie wczytanych jako łańcuchy
# (zostają tylko łańcuchy zakończone znakiem 'TER')
# -------------------------------------------------------------------------------------

            if self.data_correctnes_flag == True:
                try:
                    chains_to_remove = []
                    for index in range(len(self.CHAINS)):
                        if type(self.CHAINS[index].RESIDUES[-1]) != Termination_symbol:
                            chains_to_remove.append(index)
                        else:
                            self.CHAINS[index].RESIDUES.pop(-1)

                    chains_to_remove.reverse()
                    for index in chains_to_remove:
                        self.CHAINS.pop(index)

                except:
                    self.correctnes_flag = False
                    print(Colors.RED + f'Plik: \"{pdb_file_path}\" zawiera bledy w oznaczeniach lancuchow :O' + Colors.END)

# -------------------------------------------------------------------------------------
# Wyodrębnienie struktur drugorzędowych (wraz ze sprawdzeniem poprawności zakresów)
# -------------------------------------------------------------------------------------

            if self.data_correctnes_flag == True:
                try:
                    for chain in self.CHAINS:
                        residues_seq_numbers = [residue.residue_sequence_number for residue in chain.RESIDUES]
                        for index_h in helices_indexes_raw:
                            if chain.chain_name == index_h[1]:
                                chain.push_helix_indexes(index_h[0], index_h[2], index_h[3])
                                start_index = residues_seq_numbers.index(index_h[2])
                                end_index = residues_seq_numbers.index(index_h[3])
                                if (start_index < end_index) and (start_index < len(residues_seq_numbers)) and (end_index < len(residues_seq_numbers)):
                                    chain.push_helix(index_h[1], index_h[0], chain.RESIDUES[start_index : end_index + 1])
                                else:
                                    raise IndexError

                        for index_s in sheets_indexes_raw:
                            if chain.chain_name == index_s[1]:
                                chain.push_sheet_indexes(index_s[0], index_s[2], index_s[3])
                                start_index = residues_seq_numbers.index(index_s[2])
                                end_index = residues_seq_numbers.index(index_s[3])
                                if (start_index < end_index) and (start_index < len(residues_seq_numbers)) and (end_index < len(residues_seq_numbers)):        
                                    chain.push_sheet(index_s[1], index_s[0], chain.RESIDUES[start_index : end_index + 1])
                                else:
                                    raise IndexError

                except:
                    self.correctnes_flag = False
                    print(Colors.RED + f'Plik: \"{pdb_file_path}\" zawiera bledy w oznaczeniach alfa-helis lub beta-kartek :O' + Colors.END)

        except FileNotFoundError:
            self.correctnes_flag = False
            print(Colors.RED + f'Plik: \"{pdb_file_path}\" nie zostal znaleziony :/' + Colors.END)

# -------------------------------------------------------------------------------------
# Usunięcie częściowo wczytanych danych w przypadku wystąpienia błędu
# -------------------------------------------------------------------------------------

        if self.data_correctnes_flag == False:
            self.CHAINS.clear()

# -------------------------------------------------------------------------------------
# Funkcje pomocnicze
# -------------------------------------------------------------------------------------

    def get_residues_number(self):
        residues_number = 0
        for chain in self.CHAINS:
            residues_number += chain.get_residues_number()

        return residues_number

    def get_atoms_number(self):
        atoms_number = 0
        for chain in self.CHAINS:
            atoms_number += chain.get_atoms_number()

        return atoms_number

