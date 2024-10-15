# Contains some helper functions to format expressions 

class Str_Helper:

    def format(coefficients, variables):

        """
        
        Returns a formatted string representing the expression 
        c_0 + c_1 * x_1 + c_2 * x_2 + ... + c_n * x_n
        where (c_0, c_1, ..., c_n) is the coefficients vector and 
        (x_1, x_2, ..., x_n) is the variables vector. 

        """

        f = ""
        if coefficients[0] != 0:
            f += str(coefficients[0])

        for i in range(1, len(coefficients)):
            v = variables[i - 1]
            c = coefficients[i]

            if c == 0:
                continue
            elif c > 0:
                if len(f) != 0: f += " + "
            else:
                f += " - " if len(f) != 0 else "-"

            if c == 1 or c == -1:
                f += v
            else:
                f += str(abs(c)) + v
        return f
    