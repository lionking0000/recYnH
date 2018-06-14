import sys

class Command:
    def __init__(self):
        self.command_list = []
        self.command_dic = {}

    def AddCommand(self, name, function):
        self.command_dic[ name ] = function
        self.command_list.append( name )

    def Run(self):
        if len(sys.argv) >= 2:
            command = sys.argv[1]
        else:
            command = ""

        if command in self.command_dic:
            self.command_dic[ command ]( sys.argv )
        else:
            for command in self.command_list:
                print "python", sys.argv[0], command
