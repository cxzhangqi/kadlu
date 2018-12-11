
import matlab.engine as mat

eng = mat.start_matlab()

t = eng.gcd(100.0, 80.0, nargout=3)
print(t)

eng.quit()
