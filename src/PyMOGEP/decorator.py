#-*-coding:utf-8-*-
'''
@author: Hung-Hsin Chen
@mail: chenhh@par.cse.nsysu.edu.tw
@license: GPLv2
'''
import functools

def symbol(sym):
    '''
    Decorator that assigns a symbol to a function. 
    The symbol is stored in the function.symbol attribute.
    @param sym: symbol to a function
    '''
    def decorator(func):
        '''
        Attaches a symbol to a function as its 'symbol' attribute
        @param func: function to decorate
        '''
        func.symbol = sym
        return func
    return decorator


def cache(func):
    ''' 
    cache result of the class member method which has no argument.
    The return value is cached on self._{method}_cache where 
    {method} is the name of the method.
    usage:
        @cache
        def _get_something(self):
            ...
            return 'something'
    '''
    cache_name = '_%s_cache' %(func.func_name)

    @functools.wraps(func)
    def decorator(self):
        '''Assigns a cache attribute to self on demand'''
        try:
            return getattr(self, cache_name)
            
        except AttributeError:
            # Haven't cached anything yet
            setattr(self, cache_name, func(self))
            return getattr(self, cache_name)

    return decorator

# 
# def memory(func):
#     '''
#     cache result of the class member method which has exact one argument. 
#     self._{method}_memory where {method} is the  name of the method.
#     
#     Note that the arg must be hashable, thus lists can't be memoized.  
#     The name of the memoized attribute is stored on the method 
#     itself as func.memory.
#     usage:
#         @memoize
#         def _compute_something(self, arg):
#             ...
#             return 'something'
#     '''
#     func.memory = memory_name = '_%s_memory' %( func.func_name)
#     
#     @functools.wraps(func)
#     def decorator(self, key):
#         '''Assigns a memo hash to self on demand'''
#         try:
#             memo = getattr(self, memory_name)
#         except AttributeError:
#             # Haven't memoized anything yet
#             memo = {}
#             setattr(self, memory_name, memo)
#         
#         try:
#             return memo[key]
#         except KeyError:
#             # Haven't seen this key yet
#             memo[key] = results = func(self, key)
#             return results
#     return decorator
