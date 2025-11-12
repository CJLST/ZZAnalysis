import Mela


def initializeMELA(runMELA, year):
    """initializes a pointer to a MELA object."""
    if runMELA == True:
        sqrts = 13. 
        if year >= 2022:
            sqrts = 13.6
        import sys, os
        devnull = os.open(os.devnull, os.O_WRONLY)
        stdout_fd = sys.stdout.fileno()
        saved_stdout = os.dup(stdout_fd)
        os.dup2(devnull, stdout_fd)
        m = Mela.Mela(sqrts, 125, Mela.VerbosityLevel.ERROR) 
        m.setCandidateDecayMode(Mela.CandidateDecayMode.CandidateDecay_ZZ)  
        os.dup2(saved_stdout, stdout_fd)
        os.close(devnull)
        os.close(saved_stdout)
        print(f"***initializeMELA: created Mela({sqrts:.1f},125,TVar.CandidateDecay_ZZ)", flush=True)
    else: 
        m = None
    return m



def check_enum(entry, enum):
    found = False
    i = 0
    possible_value = tuple(enum.__members__.keys())
    mapping = enum.__members__
    while (not found) and (i < len(possible_value)):
        if entry.lower() == possible_value[i].lower():
            found = True
            entry = possible_value[i]
        i += 1
    if not found:
        possible_value = tuple(map(str.lower, possible_value))
        errortext = "Unknown matrix element given!"
        errortext += "\nThe following are valid matrix elements"
        errortext += "\n" + "\n".join(possible_value)
        errortext = print_msg_box(errortext, title="ERROR")
        raise ValueError("\n" + errortext)
    return mapping[entry]




def print_msg_box(msg, indent=1, width=0, title=""):
    """returns message-box with optional title.
    Ripped from https://stackoverflow.com/questions/39969064/how-to-print-a-message-box-in-python
    
    Parameters
    ----------
    msg : str
        The message to use
    indent : int, optional
        indent size, by default 1
    width : int, optional
        box width, by default 0
    title : str, optional
        box title, by default ""
    """
    
    lines = msg.split('\n')
    space = " " * indent
    if not width:
        width = max(map(len, lines + title.split("\n")))
    box = f'╔{"═" * (width + indent * 2)}╗\n'  # upper_border
    if title:
        box += f'║{space}{title:<{width}}{space}║\n'  # title
        box += f'║{space}{"-" * len(title):<{width}}{space}║\n'  # underscore
    box += ''.join([f'║{space}{line:<{width}}{space}║\n' for line in lines])
    box += f'╚{"═" * (width + indent * 2)}╝'  # lower_border
    return box
