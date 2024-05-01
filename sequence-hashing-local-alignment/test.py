import unittest
from main import find_matches
def test():
  unittest.main()

class LocalAlignTester(unittest.TestCase):
  def test_findsbasic(self):
    s1 = ("a" * 30) + "dcb"
    s2 = ("a" * 30) + "bcd"
    matches = find_matches(s1, s2)
    eq = [(0,0)]
    self.assertEqual(len(eq), len(matches))
    for i, m in enumerate(matches):
      self.assertEqual(eq[i][0],m[0])
      self.assertEqual(eq[i][1],m[1])

  def test_findsame(self):
    s1 = ("a" * 30) + "cda"
    s2 = ("a" * 30) + "cda"
    matches = find_matches(s1, s2)
    eq = [(x,x) for x in range(4)]
    self.assertEqual(len(eq), len(matches))
    for i, m in enumerate(matches):
      self.assertEqual(eq[i][0],m[0])
      self.assertEqual(eq[i][1],m[1])

if __name__ == '__main__':
  test()