import constants.variables as const
from utils import util


'''
Функция логирования данных в консоль
параметр - текст
value - значение
force - вывести в консоль принудительно
'''
def log(text, value="", force=False):
    if const.LOG or force:
        print(text)
        print(value)
        if util.is_not_blank(value):
            print('')