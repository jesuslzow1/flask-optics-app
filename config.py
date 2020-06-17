import os

class Config(object):
    SECRET_KEY = os.environ.get('SECRET_KEY') or 'you-will-never-an-optics-user'
    TEMPLATES_AUTO_RELOAD = True
    #DEBUG = True