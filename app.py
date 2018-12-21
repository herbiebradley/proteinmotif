from flask import Flask, render_template, url_for, request, send_file

from motif_search import match_list_from_uniprot

# Create the application object:
app = Flask(__name__)
# Use decorators to link the function to a url:
@app.route("/", methods=['GET', 'POST'])
def home():
    # Post method is called if submit button is pressed.
    if request.method == 'POST':
        # Get the contents of the "uniprot" text field:
        ac_code = request.form["uniprot"]
        # Get processed results or validation error:
        results, error = match_list_from_uniprot(ac_code)
        # Set uniprot url of protein:
        prot_url = "https://uniprot.org/uniprot/" + ac_code
    else:
        results, error, ac_code, prot_url = None, None, None, None
    # Render home.html, passing in the results, validation error, accession code, and uniprot url:
    return render_template('home.html', results=results, error=error, ac_code=ac_code, prot_url=prot_url)

@app.route('/return-file/')
def return_file():
    try:
        return send_file('results.csv', as_attachment=True)
    except Exception as e:
        return str(e)

# Create about page:
@app.route("/about")
def about():
    return render_template('about.html', title='About')

# Start the server with the 'run()' method:
if __name__ == '__main__':
    app.run(debug=True)
