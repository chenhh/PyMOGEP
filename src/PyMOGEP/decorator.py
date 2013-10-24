#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''

def symbol(symb):
    '''
    Decorator that assigns a symbol to a function. 
    The symbol is stored in the function.symbol attribute.
    @param symb: symbol to a function, typically one character
    '''
    def decorator(func):
        '''
        Attaches a symbol to a function as its 'symbol' attribute
        @param func: function to decorate
        '''
        func.symbol = symb
        return func
    return decorator

if __name__ == '__main__':
    pass