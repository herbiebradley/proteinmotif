# import the Flask class from the flask module
from flask import Flask, render_template, url_for

# create the application object
app = Flask(__name__)

# use decorators to link the function to a url
@app.route("/")
@app.route("/home")
def home():
    return render_template('home.html')

@app.route("/help")
def help():
    return render_template('help.html', title='Help')

# start the server with the 'run()' method
if __name__ == '__main__':
    app.run(debug=True)
