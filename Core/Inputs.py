#!/usr/bin/env python
# coding: utf-8
# **********************************************
# Python modules to work with 
# **********************************************

import sys
# Definitions for user inputs:::
class Inputs:
    """
    This class takes different variables to convert it in strings or integers or float depends on the needs.
    """

    def __init__(self, message):
        self.mess = message

    def input_number(self):

        while True:
            try:
                user_input = int(input(self.mess))
            except ValueError:
                print("\t> Not an integer!, please try again.")
                continue
            else:
                return user_input
                break

    def switcher(self):

        while True:
            question_switcher = input(self.mess)
            if question_switcher not in ('a', 'i', 'A', 'I'):
                print("\t> Not an appropriate choice. Please try again!")
            else:
                break
        return question_switcher
    def switcher2(self):

        while True:
            question_switcher = input(self.mess)
            if question_switcher not in ('YES', 'Y', 'yes', 'y', 'NO', 'N', 'no', 'n', "Quit", "QUIT", "QUit", "QUIt", "quit", "Q", "q", "END", "end"):
                print("\t> Not an appropriate choice. Please try again!")
            else:
                break
        return question_switcher

    def main_switcher(self):

        while True: 
            try:
                #question_switcher = int(input(self.mess))
                question_switcher = (input(self.mess))
                if question_switcher in ["Quit", "QUIT", "QUit", "QUIt", "quit", "Q", "q", "END", "end"]:
                    print("\t Leaving program...")
                    break
                else:
                    question_switcher = int(question_switcher)
                    #if question_switcher not in (1, 2, 3, 4, 5):
                    if question_switcher not in (1, 2, 3, 4, 5):
                        print("\t> Not an appropriate choice. Please try again!")
                        continue
                    else:
                        if question_switcher == 5:
                           print("\t This is a unfinished feature included in all THERM21")
                           print("\t However, not recommended to use it if you don't have")
                           print("\t any experience with python coding")
                        break
            except:
                print("\t> Not an appropriate choice. Please try again!")
                continue
        return question_switcher

    def main_switcher2(self):

        while True: 
            try:
                #question_switcher = int(input(self.mess))
                question_switcher = (input(self.mess))
                if question_switcher in ["Quit", "QUIT", "QUit", "QUIt", "quit", "Q", "q", "END", "end"]:
                    print("\t Leaving program...")
                    break
                else:
                    question_switcher = int(question_switcher)
                    #if question_switcher not in (1, 2, 3, 4, 5):
                    if question_switcher not in (0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, "all"):
                        print("\t> Not an appropriate choice. Please try again!")
                        continue
                    else:
                    
                        break
            except:
                print("\t> Not an appropriate choice. Please try again!")
                continue
        return question_switcher

    def plot_switcher(self):

        while True: 
            try:
                question_switcher = (input(self.mess))
                if question_switcher in ["Quit", "QUIT", "QUit", "QUIt", "quit", "Q", "q", "END", "end"]:
                    print("\t Leaving program...")
                    break
                else:
                    question_switcher = int(question_switcher)
                    if question_switcher not in (0, 1, 2):
                        print("\t> Not an appropriate choice. Please try again!")
                        continue
                    else:
                        break
            except:
                print("\t> Not an appropriate choice. Please try again!")
                continue
        return question_switcher

"""
    def main_switcher(self):

        while True: 
            question_switcher = int(input(self.mess))
            if question_switcher not in (1, 2, 3, 4, 5):
                print("\t> Not an appropriate choice. Please try again!")
            else:
                break
        return question_switcher

"""

