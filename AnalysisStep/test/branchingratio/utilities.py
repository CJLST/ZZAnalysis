import contextlib, os, functools

@contextlib.contextmanager
def cd(newdir):
  """http://stackoverflow.com/a/24176022/5228524"""
  prevdir = os.getcwd()
  os.chdir(os.path.expanduser(newdir))
  try:
    yield
  finally:
    os.chdir(prevdir)

class TFile(object):
  def __init__(self, *args, **kwargs):
    self.args, self.kwargs = args, kwargs
  def __enter__(self):
    import ROOT
    self.__tfile = ROOT.TFile.Open(*self.args, **self.kwargs)
    if not self.__tfile or self.__tfile.IsZombie(): return None
    return self.__tfile
  def __exit__(self, *err):
    if self.__tfile:
      self.__tfile.Close()

def cache(function):
  cache = {}
  @functools.wraps(function)
  def newfunction(*args, **kwargs):
    try:
      return cache[args, tuple(sorted(kwargs.iteritems()))]
    except TypeError:
      print(args, tuple(sorted(kwargs.iteritems())))
      raise
    except KeyError:
      cache[args, tuple(sorted(kwargs.iteritems()))] = function(*args, **kwargs)
      return newfunction(*args, **kwargs)
  return newfunction

CJLSTproduction = {
  2016: "180121",
  2017: "180416",
}
