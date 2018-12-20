from flask import Flask, render_template, url_for, request, redirect

from motif_search import match_list_from_uniprot

# create the application object
app = Flask(__name__)

# use decorators to link the function to a url
@app.route("/", methods=['GET', 'POST'])
def home():
    results = None
    if request.method == 'POST':
        input = request.form["uniprot"]
        results = match_list_from_uniprot(input)

    return render_template('home.html', results=results)

@app.route("/help")
def help():
    return render_template('help.html', title='Help')

# start the server with the 'run()' method
if __name__ == '__main__':
    app.run(debug=True)
